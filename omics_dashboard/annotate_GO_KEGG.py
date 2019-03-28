import pandas as pd
from omics_dashboard import kegg, geneontology as go
from omics_dashboard.utils import filter_de_results
import numpy as np
from functools import partial
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr


def run_GO_annotation(taxid, de_df, alpha=0.05, logFC=0.5, pval=0.05, fdr=0.05, multi_fdr=0.05):
    def GO_anno(de_):
        de_data = go.id_to_genelist(GeneID2nt, de_)
        de_data_filtered = filter_de_results(de_data, logFC, pval, fdr)
        de_data_to_GO = {'up': go.get_ID2nt(GeneID2nt, de_data_filtered, 'up'),
                         'down': go.get_ID2nt(GeneID2nt, de_data_filtered, 'down')}

        up_genes = de_data_to_GO['up']
        down_genes = de_data_to_GO['down']
        goea_results = {'up': [r for r in goeaobj.run_study(up_genes) if r.p_fdr_bh < multi_fdr],
                        'down': [r for r in goeaobj.run_study(down_genes) if r.p_fdr_bh < multi_fdr]}

        cols = ['GO', 'NS', 'p_fdr', 'count', 'terms']
        up_res = goea_results['up']
        down_res = goea_results['down']
        up_df = pd.DataFrame([[k.GO, k.NS, k.p_fdr_bh, k.study_count, k.goterm.name] for k in up_res], columns=cols)
        down_df = pd.DataFrame([[k.GO, k.NS, k.p_fdr_bh, k.study_count, k.goterm.name] for k in down_res], columns=cols)
        return up_df, down_df

    if isinstance(de_df, pd.DataFrame):
        goeaobj, GeneID2nt = go.build_network(taxid, alpha)
        up_df, down_df = GO_anno(de_df)
        goea_df = {'up': up_df,
                   'down': down_df}
        return goea_df
    else:
        goeaobj, GeneID2nt = go.build_network(taxid, alpha)
        goea_df = {}
        for study in de_df:
            up_df, down_df = GO_anno(de_df[study])
            goea_df[study] =  {'up': up_df,
                               'down': down_df}
        return goea_df


def goanna(study, taxid=511145, logFC=0.5, pval=0.05, fdr=0.05, multi_fdr=0.05):
    name = study[0]
    df = study[1]
    goeaobj, GeneID2nt = go.build_network(taxid, multi_fdr)
    de_data = go.id_to_genelist(GeneID2nt, df)
    de_data_filtered = filter_de_results(de_data, logFC, pval, fdr)
    de_data_to_GO = {'up': go.get_ID2nt(GeneID2nt, de_data_filtered, 'up'),
                     'down': go.get_ID2nt(GeneID2nt, de_data_filtered, 'down')}

    up_genes = de_data_to_GO['up']
    down_genes = de_data_to_GO['down']
    goea_results = {'up': [r for r in goeaobj.run_study(up_genes) if r.p_fdr_bh < multi_fdr],
                    'down': [r for r in goeaobj.run_study(down_genes) if r.p_fdr_bh < multi_fdr]}

    cols = ['GO', 'NS', 'p_fdr', 'count', 'terms']
    up_res = goea_results['up']
    down_res = goea_results['down']
    up_df = pd.DataFrame([[k.GO, k.NS, k.p_fdr_bh, k.study_count, k.goterm.name] for k in up_res], columns=cols)
    down_df = pd.DataFrame([[k.GO, k.NS, k.p_fdr_bh, k.study_count, k.goterm.name] for k in down_res], columns=cols)
    return ([up_df, down_df], name)


def applyParallel(func, groups, cores=None):
    from multiprocessing import Pool, cpu_count
    if cores:
        cores = cores
    else:
        cores = cpu_count()
    with Pool(cores) as p:
        ret_list = p.map(func, groups)
    return ret_list


def run_goanna(studies, taxid=511145, logFC=0.5, pval=0.05, fdr=0.05, multi_fdr=0.05, cores=None):
    func = partial(goanna, taxid=taxid, logFC=logFC, pval=pval, fdr=fdr, multi_fdr=multi_fdr)
    dfs = applyParallel(func, studies.items(), cores)
    go_dfs = {}
    for study in dfs:
        go_dfs[study[1]] = {'up':   study[0][0],
                            'down': study[0][1]}
    return go_dfs


def run_KEGG_annotation(taxid, de_df, species='eco', logFC=0.5, pval=0.05, fdr=0.05):
    symbol_to_entrez = kegg.symbol_to_entrez(taxid)

    de_data = {}
    for de_gene, row in de_df.iterrows():
        dat = row.values
        for gene in symbol_to_entrez:
            if gene.lower() == de_gene.lower():
                # gene, logFC, Pval, FDR
                de_data[symbol_to_entrez[gene]] = [dat[0], dat[3], dat[4]]

    de_data = filter_de_results(de_data, logFC, pval, fdr)

    de_kegg = {'up': {},
               'down': {}}

    importr('clusterProfiler')
    for sign in ['up', 'down']:
        if len(de_data[sign]) != 0:
            if len(de_data[sign]) == 1:
                kegg_str = "c('{0}')".format(list(de_data[sign])[0])
            else:
                kegg_str = "c{0}".format(tuple(['{0}'.format(k) for k in de_data[sign]]))

            r_cmd = 'z <- enrichKEGG(gene={0}, organism = "{1}", pvalueCutoff = {2})'.format(kegg_str, species, pval)
            r(r_cmd)

            try:
                res = r('z@result')
            except Exception as e:
                print('ERROR: enrichKEGG failed -> ', e)
                res = None

            if res:
                de_kegg[sign] = pandas2ri.rpy2py_dataframe(res)
            else:
                de_kegg[sign] = pd.DataFrame()
    return de_kegg


def kegg_anno(study, taxid=511145, species='eco', logFC=0.5, pval=0.05, fdr=0.05):
    name = study[0]
    df = study[1]

    symbol_to_entrez = kegg.symbol_to_entrez(taxid)

    de_data = {}
    for de_gene, row in df.iterrows():
        dat = row.values
        for gene in symbol_to_entrez:
            if gene.lower() == de_gene.lower():
                # gene, logFC, Pval, FDR
                de_data[symbol_to_entrez[gene]] = [dat[0], dat[3], dat[4]]

    de_data = filter_de_results(de_data, logFC, pval, fdr)

    de_kegg = {'up': {},
               'down': {}}
    importr('clusterProfiler')
    for sign in ['up', 'down']:
        if len(de_data[sign]) != 0:
            if len(de_data[sign]) == 1:
                kegg_str = "c('{0}')".format(list(de_data[sign])[0])
            else:
                kegg_str = "c{0}".format(tuple(['{0}'.format(k) for k in de_data[sign]]))

            r_cmd = 'z <- enrichKEGG(gene={0}, organism = "{1}", pvalueCutoff = {2})'.format(kegg_str, species, pval)
            r(r_cmd)

            try:
                res = r('z@result')
            except Exception as e:
                print('ERROR: enrichKEGG failed -> ', e)
                res = None

            if res:
                de_kegg[sign] = pandas2ri.rpy2py_dataframe(res)
            else:
                de_kegg[sign] = pd.DataFrame()
    return (de_kegg, name)


def run_kegg(studies, taxid=511145, species='eco', logFC=0.5, pval=0.05, fdr=0.05, cores=None):
    func = partial(kegg_anno, taxid=taxid, species=species, logFC=logFC, pval=pval, fdr=fdr)
    dfs = applyParallel(func, studies.items(), cores)
    kegg_dfs = {}
    for study in dfs:
        kegg_dfs[study[1]] = study[0]
    return kegg_dfs


def combine_annotations(go_df, kegg_df):
    all_go_terms = []
    all_kegg_terms = []
    studies = set(list(go_df) + list(kegg_df))
    for study in studies:
        if study:
            for sign in ['up', 'down']:
                if study in kegg_df:
                    if sign in kegg_df[study]:
                        if len(kegg_df[study][sign]) > 0:
                            all_kegg_terms += list(kegg_df[study][sign]['ID'].values)
                if study in go_df:
                    if sign in go_df[study]:
                        if len(go_df[study][sign]) > 0:
                            all_go_terms += list(go_df[study][sign]['GO'].values)

    total_terms = len(set(all_go_terms)) + len(set(all_kegg_terms))

    # Initialize the dataframe that will hold all annotation results
    go_anno_padj = {'up': pd.DataFrame(np.zeros((len(studies), total_terms + 1))),
                    'down': pd.DataFrame(np.zeros((len(studies), total_terms + 1)))}
    go_anno_padj['up'].columns = ['study'] + list(set(all_kegg_terms)) + list(set(all_go_terms))
    go_anno_padj['down'].columns = ['study'] + list(set(all_kegg_terms)) + list(set(all_go_terms))
    go_anno_padj['up']['study'] = kegg_df.keys()
    go_anno_padj['down']['study'] = kegg_df.keys()

    # Populate the dataframe holding all annotation results
    for study in studies:
        if study:
            for sign in ['up', 'down']:
                if study in kegg_df:
                    if sign in kegg_df[study]:
                        if not isinstance(kegg_df[study][sign], dict):
                            if 'p.adjust' in kegg_df[study][sign].columns:
                                t = go_anno_padj[sign].loc[go_anno_padj[sign]['study'] == study].merge(
                                    kegg_df[study][sign][['p.adjust']].T, how='right')
                                t['study'] = study
                                go_anno_padj[sign].loc[go_anno_padj[sign]['study'] == study] = t.values
                if study in go_df:
                    if sign in go_df[study]:
                        if not go_df[study][sign].empty:
                            go_df[study][sign].index = go_df[study][sign]['GO']
                            t = go_anno_padj[sign].loc[go_anno_padj[sign]['study'] == study].merge(
                                go_df[study][sign][['p_fdr']].T, how='right')
                            t['study'] = study
                            go_anno_padj[sign].loc[go_anno_padj[sign]['study'] == study] = t.values

    c = np.zeros((len(studies), total_terms))
    for e, i in go_anno_padj['up'].iterrows():
        c[e] = i.values[1:]
    for e, i in go_anno_padj['down'].iterrows():
        if not c[e].any():  # need to check if the row has 'up' vals
            c[e] = -i.values[1:]
        else:
            # 'up' value exists
            for n, j in enumerate(i.values[1:]):
                if 0 < j <= 0.05 and 0 < c[e][n] <= 0.05:
                    # down and up are significant
                    c[e][n] = np.inf
                elif 0 < j <= 0.05 and c[e][n] > 0.05:
                    # down is sig and up is not
                    c[e][n] = -j
                elif j > 0.05 and 0 < c[e][n] <= 0.05:
                    # down is not sig and up is
                    continue
    c[c == 0] = np.nan
    c_df = pd.DataFrame(c,
                        columns=go_anno_padj['up'].columns[1:],
                        index=go_anno_padj['up'].loc[:, 'study'])
    c_df = c_df.dropna(thresh=1, axis=0).dropna(thresh=1, axis=1)
    return c_df


def functional_annotations(go_df, kegg_df):
    go_kegg_to_func = {}
    anno = go_df
    for study in anno:
        for sign in anno[study]:
            if not isinstance(anno[study][sign], dict) and not anno[study][sign].empty:
                for i, row in anno[study][sign].iterrows():
                    term = row[0]
                    loc = row.name
                    annos = list(anno[study][sign][['terms', 'NS']].loc[loc])
                    annos = [annos[0].replace("\'", '"'), annos[1]]
                    go_kegg_to_func[term] = annos
    anno = kegg_df
    for study in anno:
        for sign in anno[study]:
            if not isinstance(anno[study][sign], dict) and not anno[study][sign].empty:
                for i, row in anno[study][sign].iterrows():
                    term = row[0]
                    loc = row.name
                    desc = anno[study][sign]['Description'].loc[loc]
                    desc = desc.replace('\'', '')
                    go_kegg_to_func[term] = [desc, 'KEGG']
    return go_kegg_to_func
