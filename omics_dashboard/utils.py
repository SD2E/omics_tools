import pandas as pd
import os
import numpy as np
from itertools import combinations
from collections import OrderedDict


def create_tempfile(df):
    name = 'edgeR_matfile.csv'
    df.to_csv(name, index=False, encoding='utf-8')
    return name


def group_by_factors(dataframe, factors):
    passed_qc_vals = dataframe.drop(factors, axis=1).values
    passed_qc_index = dataframe[factors].values
    passed_qc_cols = dataframe.drop(factors, axis=1).columns
    group_df = pd.DataFrame(passed_qc_vals, index=[list(x) for x in passed_qc_index.T], columns=passed_qc_cols)
    groups = {}
    labels = []
    group_iter = 1
    for i, (e, d) in enumerate(group_df.iterrows()):
        if e in groups:
            groups[e][1].append(i)
            groupnum = groups[e][0]
        else:
            groups[e] = [group_iter, [i]]
            groupnum = group_iter
            group_iter += 1
        labels.append(groupnum)
    return labels


def load_DE_results(de_results_dir):
    de_dfs = {}
    for f in os.listdir(de_results_dir):
        if '-vs' in f:
            fname = f.rstrip('.txt')
            try:
                de_dfs[fname] = pd.read_csv(de_results_dir + f,
                                            delimiter=' ',
                                            dtype={'logFC': 'float',
                                                   'logCPM': 'float',
                                                   'F': 'float',
                                                   'PValue': 'float',
                                                   'FDR': 'float'},
                                            float_precision='round_trip')
            except Exception as e:
                print('ERROR: ', e, f)
    return de_dfs


def prepare_dataframe(dataframe, factors=None, metadata=None, transpose=False):
    sfx = dataframe.split('.')[-1]
    if 'csv' in sfx or 'txt' in sfx:
        df = pd.read_csv(dataframe, index_col=0, header=None)
    elif 'xlsx' in sfx:
        df = pd.read_excel(dataframe)
    elif 'json' in sfx:
        df = pd.read_json(dataframe)
    else:
        raise('ERROR: Could not load dataframe ending with {0}'.format(sfx))

    if transpose:
        df = df.T

    genelist = []
    if metadata:
        if factors:
            try:
                df = df.drop(metadata, axis=1)
                genelist = list(df.drop(factors, axis=1).columns)
            except KeyError as e:
                raise('ERROR: did not find columns in dataframe, try transpose matrix.')
        else:
            try:
                df = df.drop(metadata, axis=1)
                genelist = list(df.drop(metadata, axis=1).columns)
            except KeyError as e:
                raise('ERROR: did not find columns in dataframe, try transpose matrix.')
    else:
        if factors:
            try:
                genelist = list(df.drop(factors, axis=1).columns)
            except KeyError as e:
                raise ('ERROR: did not find columns in dataframe, try transpose matrix.')

    if genelist:
        df[genelist] = df[genelist].apply(pd.to_numeric, errors='coerce')

    df = df.apply(lambda x: x.str.replace('[^A-Za-z0-9\s]+', '') if x.dtype == "object" else x)

    return df


def get_base_comparisons(dataframe, base_factor):
    return list(combinations(set(dataframe[base_factor].values), 2))


def get_factor_categories(sub_factors, dataframe):
    factor_to_categories = OrderedDict()
    sub_factors = sorted(sub_factors)
    for factor in sub_factors:
        set_factors = set(dataframe[factor])
        factor_to_categories[factor] = [str(x).strip() for x in set_factors]
    return factor_to_categories


def filter_de_results(de_data, logFC=1.0, pval=0.05, fdr=0.15):
    de_results = {'up': {}, 'down': {}}
    #logFC, Pval, FDR
    filter_de =   {x: y for x, y in de_data.items() if y[1] < pval}
    filter_de =   {x: y for x, y in filter_de.items() if y[2] < fdr}
    filter_up =   {x: y for x, y in filter_de.items() if y[0] > logFC}
    filter_down = {x: y for x, y in filter_de.items() if y[0] < -1*logFC}

    de_results['up'] = list(filter_up.keys())
    de_results['down'] = list(filter_down.keys())
    return de_results


def log10_conv(x, max_cut=9):
    x = x.fillna(value=0)
    x[x > 0] = -1*np.log10(x)
    x[x < 0] = np.log10(-1*x)
    if max_cut:
        x[x > max_cut]  = max_cut
        x[x < -max_cut] = -max_cut
    return x


def create_clustergrammer_matrix_file(annotated_dataframe, factor_to_categories, terms=False, fname='clustergrammer'):
    # Rework the first columns to allow groupby factors
    cg_groups = {}
    for studynumber, study in annotated_dataframe.iterrows():
        cg_groups[studynumber] = {}
        studienames = studynumber.split('-vs-')
        factors = studienames[0].split('_')
        for i, fact in enumerate(factors[1:]):
            category = list(factor_to_categories.items())[i]
            if fact not in category[1]:
                print('ERROR: Factor cannot associate with correct category.')
            cg_groups[studynumber][category[0]] = fact

    cg_col = []
    factor_order = sorted([cg_groups[x].keys() for x in cg_groups], key=lambda x:len(x))[-1]
    factor_order = list(factor_order)
    for x in annotated_dataframe.index:
        study = cg_groups[x]
        study_munge = {x.lstrip(): v.replace(':', '') for x, v in study.items()}
        k = str(study_munge)
        k = k.replace('{', '').replace('}', '').replace('\'', '')
        k = k.split(',')
        k = [l.lstrip() for l in k]

        a_f = [i for i in factor_order if i not in study.keys()]
        for i in a_f:
            k.append('{0}: Null'.format(i.strip()))

        kk = []
        for z in factor_order:
            for zz in k:
                zzz = zz.split(':')[0]
                if z in zzz:
                    kk.append(zz)
        kk = ', '.join(['\'{0}\''.format(h) for h in kk])
        col = '(\'{0}\', {1})'.format(x.replace(':', ''), kk)
        cg_col.append(col)
    annotated_dataframe = annotated_dataframe.set_index([cg_col], append=True).reset_index(level='study', drop=True)

    # Rework the first rows to allow groupby annotations
    if terms:
        anno_col = ['(\'Desc:{0}\', \'NS: {1}\')'.format(terms[i][0], terms[i][1]) for i in annotated_dataframe.columns]
        annotated_dataframe.columns = anno_col

    annotated_dataframe = annotated_dataframe.fillna(value=0).copy()
    annotated_dataframe = annotated_dataframe.replace(np.inf, 0).copy()
    annotated_dataframe.to_csv(fname + '.txt', sep='\t')

def filter_df(filter_string, df, pthresh=1):
    filter_ = []
    for x in df.index:
        if filter_string in x[0]:
            continue
        elif 'Null' in filter_string:
            for i in range(1,len(x)):
                for k in x[i]:
                    if 'Null' in k:
                        continue
        else:
            filter_.append(x)
    if 'Null' in filter_string:
        for x in df.index:
            for i in range(len(x)):
                if 'Null' in x[i]:
                    filter_.append(x)

    df = df.loc[filter_]

    df[(0 < df) & (df < pthresh)] = 0
    df[(0 > df) & (df > -pthresh)] = 0

    df[df==0] = np.nan
    df = df.dropna(thresh=1, axis=1).dropna(thresh=1, axis=0)
    df.fillna(value=0, inplace=True)
    return df


def subset_df(df, strains=None, pivots=None, constants=None, pthresh=1, max_p=9):
    filter_ = []
    if strains:
        for strain in strains:
            for x in df.index:
                if strain in x[0]:
                    filter_.append(x)
        df = df.loc[filter_]

    filter_ = []
    sub_factor_order = df.index[0][1:]
    if pivots:
        for pivot in pivots:
            for x in df.index:
                comp = x[0].split('-vs-')
                s1, s2 = comp[0].split('_'), comp[1].split('_')
                pivot_index = [i for i, el in enumerate(sub_factor_order) if pivot.lower() in sub_factor_order[i].lower()]

                same = [s1[i] == s2[i] for i in range(1, len(s1))]
                if not pivot_index:
                    print('ERROR: could not find category to pivot on')
                if same[pivot_index[0]] == False:
                    filter_.append(x)
        df = df.loc[filter_]

    constant_filter = []
    if constants:
        if not isinstance(constants, dict):
            print('ERROR constants must be a dictionary {\'Time\': \'8\'')

    if constants:
        for x in df.index:
            comp = x[0].split('-vs-')
            s1, s2 = comp[0].split('_'), comp[1].split('_')
            flags = len(constants) * [False]
            for z, const in enumerate(constants):
                const_index = [i for i, el in enumerate(sub_factor_order) if const.lower() in x[1 + i].lower()]
                if (s1[const_index[0] + 1].lower() == constants[const].lower() and s2[const_index[0] + 1].lower() ==
                        constants[const].lower()):
                    flags[z] = True
            if all(flags):
                constant_filter.append(x)

    if constant_filter:
        df = df.loc[constant_filter]

    df[(0 < df) & (df < pthresh)] = 0
    df[(0 > df) & (df > -pthresh)] = 0

    df[df > max_p] = max_p
    df[df < -max_p] = -max_p

    df[df == 0] = np.nan
    df = df.dropna(thresh=1, axis=1).dropna(thresh=1, axis=0)
    if df.empty:
        print('Dataframe is empty, filters too stringent.')
    df.fillna(value=0, inplace=True)
    return df
