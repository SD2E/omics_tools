# !/bin/bash
# -t 1-112
# -tc 36
# -l h_rt=3600
# -l h_vmem=4G
# -l m_mem_free=4G
# -o /btl/foundry/users/alex/20190228_novel_chassis/grid_out/
# -e /btl/foundry/users/alex/20190228_novel_chassis/grid_err/
source /broad/software/scripts/useuse
reuse UGER
infile=$(awk "NR==$SGE_TASK_ID" /btl/foundry/users/alex/20190228_novel_chassis/run_DGE/edgeR_files.txt)
sh $infile