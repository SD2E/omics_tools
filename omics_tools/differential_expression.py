from omics_tools import utils, comparison_generator
import os
import pandas as pd
from rpy2.robjects import r


def make_contrast_strings(comparison_indices, groups_array):
    contrast_strings = []
    for comp in comparison_indices:
        unrolled = comparison_indices[comp]
        indices = list(unrolled[0]) + list(unrolled[1])
        group_a = [groups_array[int(j)] for j in list(unrolled[0])]
        group_b = [groups_array[int(j)] for j in list(unrolled[1])]
        unrolled_digits = [0] * len(set(groups_array))
        for a in group_a:
            unrolled_digits[a - 1] = 1
        for b in group_b:
            unrolled_digits[b - 1] = -1
        c = 'c({})'.format(','.join(map(str, unrolled_digits)))
        contrast_strings.append((comp, c, indices))
    return contrast_strings


def format_groups_array(groups_array):
    if len(groups_array) > 300:
        formatted_cmd = 'group <- c{}'.format(tuple(groups_array[:300]))

        for i in range(300, len(groups_array), 300):
            chunk = str(tuple(groups_array[i:i + 300])).replace('(', '').replace(')', '')
            formatted_cmd += '\ngroup <- c(group, {})'.format(chunk)
        formatted_cmd += '\ngroup <- factor(group)'
    else:
        formatted_cmd = '\ngroup <- factor(c{})'.format(tuple(groups_array))
    return formatted_cmd


def make_hpc_de_files(dataframe=None, base_comparisons=None, data_frame_path=None, base_factor=['strain'],
         sub_factors=None, freedom=1, metadata=None, transpose=False, run_dir=None):

    sub_factors = sorted(sub_factors)
    if not isinstance(dataframe, pd.DataFrame):
        dataframe = utils.prepare_dataframe(data_frame_path, base_factor + sub_factors, metadata, transpose)
    dataframe.reset_index(drop=True, inplace=True)

    df_file = os.path.basename(utils.create_tempfile(dataframe))
    if not os.path.exists('./scripts/'):
        os.mkdir('./scripts/')

    if not base_comparisons:
        base_comparisons = utils.get_base_comparisons(dataframe, base_factor)
    comparison_indices = comparison_generator.generate_comparisons(dataframe, base_comparisons, base_factor,
                                                                   sub_factors, freedom)
    groups_array = utils.group_by_factors(dataframe, base_factor + sub_factors)
    contrast_strings = make_contrast_strings(comparison_indices, groups_array)

    files = open('./edgeR_files.txt', 'w')
    for e, c in enumerate(contrast_strings):
        filt0 = 'keep <- rowSums(cpm(y[, c{0}]) >1) >= {1}'.format(tuple(c[2]), len(c[2]))
        filt1 = 'y <- y[keep, ,]'
        filt2 = 'y <- estimateDisp(y, design)'
        str0 = 'fit <- glmQLFit(y, design)'
        str1 = 'qlf <- glmQLFTest(fit, contrast={})'.format(c[1])
        str2 = 'tab <- topTags(qlf, n=Inf)'
        c_ = '-vs-'.join(['_'.join(map(str, x)) for x in c[0]])
        c_ = c_.replace(' ', '')
        str3 = 'write.table(tab, file="{}.txt")'.format(c_)
        R_string = '\n'.join([filt0, filt1, filt2, str0, str1, str2, str3])
        fn = './scripts/dge_{}.r'.format(str(e))
        if run_dir:
            with open(fn, 'w') as of:
                of.write(Rscript(run_dir+df_file, groups_array, factors=base_factor + sub_factors))
                of.write(R_string)
        else:
            with open(fn, 'w') as of:
                of.write(Rscript(df_file, groups_array, factors=base_factor+sub_factors))
                of.write(R_string)
        sh_fn = './scripts/dge_{}.sh'.format(str(e))
        if run_dir:
            sh_fn = run_dir + '/scripts/dge_{}.sh'.format(str(e))
        files.write(sh_fn + '\n')
        if run_dir:
            exec_script = Exec_script(e, run_dir)
        else:
            exec_script = Exec_script(e)
        script_file = './scripts/dge_{}.sh'.format(str(e))
        with open(script_file, 'w') as of:
            of.write(exec_script)
    files.close()

    with open('./dge_array.sh', 'w') as f:
        cmd = """# !/bin/bash
# -t 1-{0}
# -tc 36
# -l h_rt=3600
# -l h_vmem=4G
# -l m_mem_free=4G
# -o /btl/foundry/users/alex/20190228_novel_chassis/grid_out/
# -e /btl/foundry/users/alex/20190228_novel_chassis/grid_err/
source /broad/software/scripts/useuse
reuse UGER
infile=$(awk "NR==$SGE_TASK_ID" /btl/foundry/users/alex/20190228_novel_chassis/run_DGE/edgeR_files.txt)
sh $infile""".format(str(len(contrast_strings)))
        f.write(cmd)
    return 0


def make_DE_cmds(dataframe=None, base_comparisons=None, base_factor=['strain'],
         sub_factors=None, freedom=1, run_dir=None, tofile=False):

    if not base_comparisons:
        base_comparisons = utils.get_base_comparisons(dataframe, base_factor)

    sub_factors = sorted(sub_factors)
    groups_array = utils.group_by_factors(dataframe, base_factor + sub_factors)
    comparison_indices = comparison_generator.generate_comparisons(dataframe, base_comparisons, base_factor,
                                                                   sub_factors, freedom)
    df_file = os.path.basename(utils.create_tempfile(dataframe))
    if run_dir:
        df_file = run_dir + df_file

    contrast_strings = make_contrast_strings(comparison_indices, groups_array)

    R_cmds = []
    for e, c in enumerate(contrast_strings):
        c_ = '-vs-'.join(['_'.join(map(str, x)) for x in c[0]])
        c_ = c_.replace(' ', '')
        name0  = '#' + c_
        edger0 = Rscript2(df_file, base_factor + sub_factors)
        edger1 = 'group <- factor(c{0})'.format(format_groups_array(groups_array))
        edger2 = 'y <- DGEList(counts=t_cts, group=group)'
        edger3 = 'y <- calcNormFactors(y)'
        edger4 = 'design <- model.matrix(~0 + group)'

        filt0  = 'keep <- rowSums(cpm(y[, c{0}]) >1) >= {1}'.format(tuple(c[2]), len(c[2]))
        filt1  = 'y <- y[keep, ,]'
        filt2  = 'y <- estimateDisp(y, design)'
        str0   = 'fit <- glmQLFit(y, design)'
        str1   = 'qlf <- glmQLFTest(fit, contrast={})'.format(c[1])
        str2   = 'tab <- topTags(qlf, n=Inf)'

        R_string = '\n'.join([name0, edger0, edger1, edger2, edger3, edger4,
                              filt0, filt1, filt2, str0, str1, str2])
        R_cmds.append(R_string)
    return R_cmds


def Exec_script(e, run_dir='.'):
    return '''#!/usr/bin/env bash
source /broad/software/scripts/useuse
use R-3.5

WD="{0}"
RESULTS=$JOB_ID
RUNDIR="$WD$RESULTS"
mkdir $RUNDIR
cd $RUNDIR

Rscript $WD/scripts/dge_{1}.r > $WD/$RESULTS/dge_{1}_log.txt
'''.format(run_dir, str(e))


def Rscript(count_fname, groups_array, factors=None):
    cmd = """suppressMessages(library(edgeR))
options(scipen=999)

x <- read.delim('{0}', sep=',', stringsAsFactors=TRUE)
""".format(count_fname)
    for factor in factors:
        cmd += """
x${0} <- factor(x${0})
""".format(factor)
    cmd += """
drops <- c{0}
counts <- x[, !(names(x) %in% drops)]
t_cts <- t(counts)
t_cts[is.na(t_cts)] <- 0

#colnames(t_cts) <- x$filename

{1}

y <- DGEList(counts=t_cts, group=group)
y <- calcNormFactors(y)
design <- model.matrix(~0 + group)
""".format(tuple(factors), format_groups_array(groups_array))
    return cmd


def Rscript2(count_fname, factors):
    cmd = """
suppressMessages(library(edgeR))
options(scipen=999)
x <- read.delim('{0}', sep=',', stringsAsFactors=TRUE)
""".format(count_fname)
    for factor in factors:
        cmd += """
x${0} <- factor(x${0})
""".format(factor)
    cmd += """
drops <- c{0}
counts <- x[, !(names(x) %in% drops)]
t_cts <- t(counts)
t_cts[is.na(t_cts)] <- 0
#colnames(t_cts) <- x$filename
""".format(tuple(factors))
    return cmd


def edger(rcmd):
    from rpy2.robjects import r
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()

    r(rcmd)
    df = r('as.data.frame(tab)')

    name = rcmd.split('\n')[0][1:]
    return (df, name)


def applyParallel(groups, func, cores=None):
    from multiprocessing import Pool, cpu_count
    if cores:
        cores = cores
    else:
        cores = cpu_count()
    p = Pool(cores)
    ret_list = []
    for i in groups:
        p.apply_async(func, args=(i,), callback=ret_list.append)
    p.close()
    p.join()
    return ret_list


def check_installs():
    install_cmd = """if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("edgeR", quietly=TRUE))
    BiocManager::install("edgeR")"""
    r(install_cmd)


def run_edgeR(rcmds, cores=None):
    check_installs()
    if isinstance(rcmds, list):
        print('Running {} pair-wise differential tests through edgeR'.format(len(rcmds)))
        dfs = applyParallel(rcmds, edger, cores)
    else:
        print('Running a pair-wise differential test through edgeR')
        dfs = [edger(rcmds)]
    de_dfs = {}
    for study in dfs:
        de_dfs[study[1]] = study[0]
    return de_dfs

