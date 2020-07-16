import pandas as pd
from clustergrammer_widget import *
import pickle
from omics_tools import differential_expression, annotate_GO_KEGG, utils, comparison_generator

ginkgo_data = '/Users/meslami/Documents/GitRepos/omics_tools/tests/Reordered_ReadCountMatrix_preCAD.csv'

ginkgo_df = pd.read_csv(ginkgo_data, index_col=0, header=None, low_memory=False).T
ginkgo_df['strain'] = [x.split('/')[-2].replace('_','') for x in ginkgo_df['strain'].values]
ginkgo_df.rename(columns={'timepoint':'Timepoint', 'temp': 'Temp'}, inplace=True)
ginkgo_df['IPTG'] = ginkgo_df['IPTG'].replace(' NA', 0)
ginkgo_df['IPTG'] = ginkgo_df['IPTG'].astype(bool)
ginkgo_df['Arabinose'] = ginkgo_df['Arabinose'].replace(' NA', 0)
ginkgo_df['Arabinose'] = ginkgo_df['Arabinose'].astype(bool)
ginkgo_df['Timepoint'] = ginkgo_df['Timepoint'].replace(regex=':hour', value='')
ginkgo_df['Temp'] = ginkgo_df['Temp'].replace(regex=':celsius', value='')
ginkgo_df.drop(['GinkgoID','R1','R2','filename','flags','library','gene_id','replicate'], axis=1, inplace=True)
for x in ['Temp', 'Timepoint']:
    ginkgo_df[x] = pd.to_numeric(ginkgo_df[x]).astype('int64')


sub_factors = ['Timepoint', 'Temp', 'Arabinose', 'IPTG']
# sub_factors = ['Timepoint', 'Temp']


DE_tests = [
    ['Strain1MG1655WT', 'Strain1MG1655WT']
]
run_dir='/Users/meslami/Desktop/test_remove'

differential_expression.make_hpc_de_files(
    dataframe=ginkgo_df,
    base_comparisons=DE_tests,
    sub_factors=sub_factors,
    run_dir=run_dir,filter_unused_base_factors=True,export_tagwise_noise=False,aggregation_flag=True)


# utils.aggregate_dataframes(run_dir,sub_factors)