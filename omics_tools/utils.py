import pandas as pd
import os
import numpy as np
from itertools import combinations, chain
from collections import OrderedDict
import json

def create_tempfile(df, fn=None):
    if fn:
        name = fn
    else:
        name = 'edgeR_matfile.csv'
    df.to_csv(name, index=False, encoding='utf-8')
    return name


def group_by_factors(dataframe, factors):
    '''
    returns a 1-based index of the R dataframe where each base factor + subfactor combination gets a unique index
    and that index is repeated across replicates which will be grouped for the analysis
    :param dataframe: counts dataframe
    :param factors: list that combines the base factor and the subfactors
    :return:
    '''
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

def remove_non_base_samples(dataframe,factors,base_factor=['strain']):
    merged = set(chain(*factors))
    new_df = dataframe[dataframe[base_factor[0]].isin(merged)]
    return new_df

def prepare_dataframe(dataframe, factors=None, metadata=None, transpose=False):
    if isinstance(dataframe,str):
        sfx = dataframe.split('.')[-1]
        if 'csv' in sfx or 'txt' in sfx:
            df = pd.read_csv(dataframe, index_col=0, header=None)
        elif 'xlsx' in sfx:
            df = pd.read_excel(dataframe)
        elif 'json' in sfx:
            df = pd.read_json(dataframe)
        else:
            raise('ERROR: Could not load dataframe ending with {0}'.format(sfx))
    df = dataframe.copy()
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

    # df = df.apply(lambda x: x.str.replace('[^A-Za-z0-9\s]+', '') if x.dtype == "object" else x)

    return df,genelist


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
    index_reordered = []
    values_reordered = []
    for e, (i, x) in enumerate(annotated_dataframe.iterrows()):
        comp = x.name.split('-vs-')
        c1 = [z for z in comp[0].split('_')]
        c2 = [z for z in comp[1].split('_')]
        cols = list(map(str, list(range(len(c1)))))

        reor_df = pd.DataFrame([c1, c2], columns=cols)
        reor_df[cols] = reor_df[cols].apply(pd.to_numeric, errors='ignore')
        reor_df.sort_values(cols, axis=0, inplace=True)

        new_comp = []
        for i, z in reor_df.iterrows():
            new_comp.append('_'.join([str(y) for y in z.values]))

        if comp[0] != new_comp[0]:
            val_idxs = np.nonzero(x.values)
            x.values[val_idxs] = -1 * x.values[val_idxs]

        new_comp = '-vs-'.join(new_comp)

        values_reordered.append(x.values)
        index_reordered.append(new_comp)
    annotated_dataframe = pd.DataFrame(values_reordered, index=index_reordered, columns=annotated_dataframe.columns)
    annotated_dataframe.index.name = 'study'

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


def subset_df(df, strains=None, pivots=None, constants=None, pthresh=1, max_p=15, other_strains=True, same=True):
    filter_ = []
    if strains:
        for x in df.index:
            z = x[0].split('-vs-')
            for strain in strains:
                if strain in z[0] or strain in z[1]:
                    filter_.append(x)

        remove_strains = []
        if not other_strains:
            for study in filter_:
                studynames = study[0].split('-vs-')
                presence = []
                for strain in strains:
                    if strain in studynames[0]:
                        presence.append(1)
                    if strain in studynames[1]:
                        presence.append(1)
                if sum(presence) != 2:
                    remove_strains.append(study)
        remove_same_strains = []
        if not same:
            for study in filter_:
                studynames = study[0].split('-vs-')
                for strain in strains:
                    if strain in studynames[0] and strain in studynames[1]:
                        remove_same_strains.append(study)

        filter_ = set(filter_) - set(remove_strains) - set(remove_same_strains)
        filter_ = list(set(filter_))
        df = df.loc[filter_]

    filter_ = []
    sub_factor_order = df.index[0][1:]
    if pivots:
        index_reordered = []
        values_reordered = []
        for e, (i, x) in enumerate(df.iterrows()):
            comp = x.name[0].split('-vs-')
            c1 = [z for z in comp[0].split('_')]
            c2 = [z for z in comp[1].split('_')]
            cols = ['study'] + [k.split(':')[0] for k in x.name[1:]]

            reor_df = pd.DataFrame([c1, c2], columns=cols)
            reor_df[cols] = reor_df[cols].apply(pd.to_numeric, errors='ignore')
            reor_df.sort_values(pivots, axis=0, inplace=True)

            new_comp = []
            for i, z in reor_df.iterrows():
                new_comp.append('_'.join([str(y) for y in z.values]))

            if comp[0] != new_comp[0]:
                val_idxs = np.nonzero(x.values)
                x.values[val_idxs] = -1 * x.values[val_idxs]

            new_comp = '-vs-'.join(new_comp)
            x_name = [new_comp] + list(x.name)[1:]

            values_reordered.append(x.values)
            index_reordered.append(tuple(x_name))
        df = pd.DataFrame(values_reordered, index=index_reordered, columns=df.columns)

        for pivot in pivots:
            for x in df.index:
                comp = x[0].split('-vs-')
                s1, s2 = comp[0].split('_'), comp[1].split('_')
                pivot_index = [i for i, el in enumerate(sub_factor_order) if
                               pivot.lower() in sub_factor_order[i].lower()]

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

def aggregate_dataframes(run_dir,subfactors,cols_to_keep=['logFC','FDR'],suffix_col='strain_2',melt=False):
    if not os.path.exists(os.path.join(run_dir,'tmp/aggregate_dict.json')):
        raise ValueError("Aggregate dict does not exist. Need to run comparison_generator.py first with aggregate_df=True")
    with open(os.path.join(run_dir,'tmp/aggregate_dict.json')) as json_file:
        agg_dict = json.load(json_file)

    genes = pd.read_csv(os.path.join(run_dir,'genes.txt'))
    genes['tmp']=1
    genes.set_index('Gene',inplace=True)
    de_dfs = {}
    noise_dfs = {}
    de_results_dir=os.path.join(run_dir,'results')
    noise_files_found=False
    for f in os.listdir(de_results_dir):
        if ('-vs' in f) and ('dispersion' not in f):
            # try:
            fname = os.path.join(de_results_dir, f)
            print('Reading:',fname)
            df = pd.read_csv(fname,
            delimiter=' ',
            dtype={'logFC': 'float',
                   'logCPM': 'float',
                   'F': 'float',
                   'PValue': 'float',
                   'FDR': 'float'},
            float_precision='round_trip')

            df = df[cols_to_keep]
            df['nlogFDR']=-1*np.log10(df['FDR'])

            df['flag_edgeRremoved']=0
            df = df.join(genes,how='outer')
            df['flag_edgeRremoved']=df['flag_edgeRremoved'].fillna(1).astype(int)
            df.fillna(0,inplace=True)
            df.drop('tmp',inplace=True,axis=1)
            df.columns = [col + '_' + agg_dict[f][suffix_col] for col in df.columns]
            group_dict = agg_dict[f].copy()
            str_dict = json.dumps(_get_dict_subfactor_overlap(group_dict,subfactors))
            print(f)
            print(str_dict)
            if str_dict not in de_dfs:
                de_dfs[str_dict] = {}
                de_dfs[str_dict]['dfs'] = [df]
                de_dfs[str_dict]['fname'] = [f]
            else:
                de_dfs[str_dict]['dfs'].append(df)
                de_dfs[str_dict]['fname'].append(f)


            # except Exception as e:
            #     print('ERROR: ', e, f)
        if ('dispersion' in f):
            f_n = f.replace('_dispersion','')
            group_dict_noise = agg_dict[f_n].copy()
            str_dict_noise = json.dumps(_get_dict_subfactor_overlap(group_dict_noise,subfactors))
            noise_files_found = True
            df_noise = pd.read_csv(os.path.join(de_results_dir, f),delimiter=' ',
                                   dtype={'genes':'str','tag_noise':'float'},
                                   float_precision='round_trip')
            df_noise.set_index('genes',inplace=True)
            df_noise.columns = [col + '_' + agg_dict[f_n][suffix_col] for col in df_noise.columns]
            if 'NAND' in f:
                print(df_noise)
            if str_dict_noise not in noise_dfs:
                noise_dfs[str_dict_noise]={}
                noise_dfs[str_dict_noise]['dfs_noise']=[df_noise]
                noise_dfs[str_dict_noise]['fname_noise']=[f]
            else:
                noise_dfs[str_dict_noise]['dfs_noise'].append(df_noise)
                noise_dfs[str_dict_noise]['fname_noise'].append(f)

    #Combine the dataframe lists
    for key in de_dfs:

        df_all = pd.concat(de_dfs[key]['dfs'], axis=1, join='outer',sort=True).fillna(0)
        json_acceptable_string = key.replace("'", "\"")
        new_dict = json.loads(json_acceptable_string)
        for key2 in new_dict:
            df_all[key2]=new_dict[key2]
        de_dfs[key]['df_all']= df_all

    massive_df = pd.concat([de_dfs[key]['df_all'] for key in de_dfs.keys()],sort=True)
    massive_df.to_csv(os.path.join(run_dir,'results/massive_df.csv'))
    #If noise files are present then output those as well:
    if noise_files_found:
        # Combine the dataframe lists
        for key in noise_dfs:
            df_noise_all = pd.concat(noise_dfs[key]['dfs_noise'], axis=1, join='outer', sort=True).fillna(0)
            json_acceptable_string = key.replace("'", "\"")
            new_dict = json.loads(json_acceptable_string)
            for key2 in new_dict:
                df_noise_all[key2] = new_dict[key2]
            noise_dfs[key]['df_noise_all'] = df_noise_all

        massive_noise_df = pd.concat([noise_dfs[key]['df_noise_all'] for key in noise_dfs.keys()], sort=True)
        massive_noise_df.to_csv(os.path.join(run_dir,'results/massive_noise_df.csv'))
    # os.remove(run_dir+'tmp/aggregate_dict.json')
    # os.rmdir(run_dir+'./tmp')

    if melt:
        pass

    return massive_df

def _get_dict_subfactor_overlap(group_dict,subfactors):
    clean_dict = {}
    for key in group_dict.keys():
        key_temp = key[:-2]
        if key_temp in subfactors:
            clean_dict[key_temp] = group_dict[key]
    return clean_dict