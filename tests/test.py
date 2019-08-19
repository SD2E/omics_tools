import pandas as pd
from clustergrammer_widget import *
import pickle
from omics_tools import differential_expression, annotate_GO_KEGG, utils, comparison_generator

ginkgo_data = '/tests/Reordered_ReadCountMatrix_preCAD.csv'

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

sub_factors = ['Timepoint', 'Temp', 'Arabinose', 'IPTG']
DE_tests = [
    ['Strain1MG1655WT', 'Strain1MG1655WT'],
    ['Strain1MG1655WT', 'Strain2MG1655GenomicPhlFGate'], #WT vs all
    ['Strain1MG1655WT', 'Strain3MG1655GenomicIcaRGate'],
    ['Strain1MG1655WT', 'Strain4MG1655GenomicNANDCircuit'],
    ['Strain4MG1655GenomicNANDCircuit', 'Strain3MG1655GenomicIcaRGate'], #NAND vs IcaR genome
    ['Strain4MG1655GenomicNANDCircuit', 'Strain2MG1655GenomicPhlFGate'], #NAND vs PhlF genome
    ['Strain4MG1655GenomicNANDCircuit', 'Strain4MG1655GenomicNANDCircuit']  #NAND vs NAND
]

groups_array = utils.group_by_factors(ginkgo_df, ['strain']+sub_factors)
comparison_indices = comparison_generator.generate_comparisons(ginkgo_df, DE_tests, ['strain'], sub_factors, 1)
contrast_strings = differential_expression.make_contrast_strings(comparison_indices, groups_array)

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
with open('go_results.pkl', 'rb') as f:
    GO_df = pickle.load(f)
anno_df = annotate_GO_KEGG.combine_annotations(GO_df, KEGG_df)
cg_file = 'new_results'
factor_to_categories = utils.get_factor_categories(sub_factors, ginkgo_df)
terms = annotate_GO_KEGG.functional_annotations(GO_df, KEGG_df)
utils.create_clustergrammer_matrix_file(
    annotated_dataframe = anno_df,
    factor_to_categories = factor_to_categories,
    terms = terms,
    fname = cg_file)
net = Network(clustergrammer_widget)
net.load_file(cg_file + '.txt')
cgdf = net.export_df()
cgdf = utils.log10_conv(cgdf)
net.load_df(cgdf)
net.cluster(dist_type='euclidean')
net.widget()