import pandas as pd
import sys
import pickle
from clustergrammer_widget import *
from omics_tools import differential_expression, annotate_GO_KEGG, utils, comparison_generator

def main(dataframe):
    ginkgo_data = dataframe

    ginkgo_df = pd.read_csv(ginkgo_data, index_col=0, header=None, low_memory=False).T
    ginkgo_df['strain'] = [x.split('/')[-2].replace('_','') for x in ginkgo_df['strain'].values]
    ginkgo_df.rename(columns={'timepoint':'Timepoint', 'temp': 'Temp'}, inplace=True)
    ginkgo_df['IPTG'] = ginkgo_df['IPTG'].replace(' NA', 0)
    ginkgo_df['IPTG'] = ginkgo_df['IPTG'].astype(bool)
    ginkgo_df['Arabinose'] = ginkgo_df['Arabinose'].replace(' NA', 0)
    ginkgo_df['Arabinose'] = ginkgo_df['Arabinose'].astype(bool)
    ginkgo_df['Timepoint'] = ginkgo_df['Timepoint'].replace(regex=':hour', value='')
    ginkgo_df['Temp'] = ginkgo_df['Temp'].replace(regex=':celsius', value='')
    ginkgo_df.drop(['GinkgoID','R1','R2','filename','flags','library','gene_id'], axis=1, inplace=True)
    for x in ['Temp', 'Timepoint']:
        ginkgo_df[x] = pd.to_numeric(ginkgo_df[x]).astype('int64')

    ginkgo_df['strain_inducer']=ginkgo_df['strain']+'_'+ginkgo_df['Arabinose'].map(str)+'_'+ginkgo_df['IPTG'].map(str)

    # sub_factors = ['Timepoint', 'Temp', 'Arabinose', 'IPTG']
    sub_factors = ['Timepoint', 'Temp']


    DE_tests = [
        ['Strain1MG1655WT_False_False', 'Strain1MG1655WT_False_True'],
        ['Strain1MG1655WT_False_False', 'Strain1MG1655WT_True_False'],
        ['Strain1MG1655WT_False_False', 'Strain1MG1655WT_True_True'],
    ]

    ginkgo_df.reset_index(inplace=True)
    groups_array = utils.group_by_factors(ginkgo_df, ['strain_inducer']+sub_factors)

    # df, base_comparisons=None, base_factor=['strain'], sub_factors=None, freedom=1,aggregation_flag=False,run_dir=None
    comparison_indices = \
        comparison_generator.generate_comparisons(ginkgo_df,
                                                  base_comparisons= DE_tests,
                                                  base_factor= ['strain_inducer'],
                                                  sub_factors=sub_factors,
                                                  freedom= 1,
                                                  aggregation_flag=True,
                                                  run_dir='.')
    contrast_strings = \
        differential_expression.make_contrast_strings(comparison_indices,
                                                      groups_array)


    r_cmds = differential_expression.make_DE_cmds(
                dataframe = ginkgo_df,
                base_comparisons = DE_tests,
                sub_factors = sub_factors)
    print('Created {0} differential tests'.format(len(r_cmds)))

    deg_results = differential_expression.run_edgeR(r_cmds, cores=24)

    KEGG_df = annotate_GO_KEGG.run_kegg(deg_results)
    with open('kegg_results.pkl', 'wb') as f:
        pickle.dump(KEGG_df, f)

    with open('kegg_results.pkl', 'rb') as f:
        KEGG_df = pickle.load(f)

    GO_df = annotate_GO_KEGG.run_go(deg_results)
    with open('go_results.pkl', 'wb') as f:
        pickle.dump(GO_df, f)

    differential_expression.make_hpc_de_files(
        dataframe=ginkgo_df,
        base_comparisons=DE_tests,
        sub_factors=sub_factors,
        run_dir='hpc_files',
        filter_unused_base_factors=True,
        export_tagwise_noise=False)

if __name__ == '__main__':
    main(sys.argv[1])
