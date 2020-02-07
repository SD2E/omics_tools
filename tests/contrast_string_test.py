import pandas as pd
import os

full_path = '/Users/meslami/Downloads/Plan-Requirements-B-subtilis-inducer-RNASeq_omics_tools.csv'
df_counts = pd.read_csv(full_path)

df_counts.drop(df_counts.columns[0],axis=1,inplace=True)
df_counts.rename({'Cuminic_acid':'ca_concentration','Vanillic_acid':'va_concentration','IPTG':'iptg_concentration','Xylose':'xylose_concentration'},inplace=True,axis=1)
for col in df_counts.columns:
    if 'concentration' in col:
        df_counts[col] = (df_counts[col].apply(float) > 0 )
        df_counts[col].astype(bool, inplace=True)

ginkgo_df = df_counts.copy()

ginkgo_df.columns = map(str.lower, ginkgo_df.columns)
ginkgo_df['strain_inducer']='wt'

for col in ginkgo_df.columns:
    if '_concentration' in col:
        inducer = col.replace('_concentration','')
        ginkgo_df['strain_inducer']=ginkgo_df['strain_inducer']+'__'+inducer+'_'+ginkgo_df[col].map(str)

factors_to_drop=['unit','arabinose','qc','strain_input_state','replicate','temperature','concentration']
for factor in factors_to_drop:
    for col in ginkgo_df.columns:
        if factor in str.lower(col):
            ginkgo_df.drop(col, axis=1, inplace=True)

sub_factors = ['timepoint']
run_dir='/Users/meslami/Desktop/DARPA/SD2/NovelChassis/Bacillus/Inducer_Prediction'
if not os.path.exists(run_dir):
    os.makedirs(run_dir)

DE_tests = [
    ['wt__ca_False__iptg_False__va_False__xylose_False', 'wt__ca_False__iptg_False__va_True__xylose_False'],
    ['wt__ca_False__iptg_False__va_False__xylose_False', 'wt__ca_False__iptg_True__va_True__xylose_False'],
    ['wt__ca_False__iptg_False__va_False__xylose_False','wt__ca_True__iptg_False__va_False__xylose_False'],
    ['wt__ca_False__iptg_False__va_False__xylose_False','wt__ca_True__iptg_True__va_False__xylose_False'],
    ['wt__ca_False__iptg_False__va_False__xylose_False', 'wt__ca_False__iptg_True__va_False__xylose_False'],
    ['wt__ca_False__iptg_False__va_False__xylose_False', 'wt__ca_False__iptg_False__va_False__xylose_True'],
]

ginkgo_df.reset_index(inplace=True,drop=True)

from omics_tools import differential_expression,  utils, comparison_generator

groups_array = utils.group_by_factors(ginkgo_df, ['strain_inducer']+sub_factors)
comparison_indices = comparison_generator.generate_comparisons(ginkgo_df,base_comparisons= DE_tests,base_factor= ['strain_inducer'], sub_factors=sub_factors,freedom= 1,aggregation_flag=True,run_dir=run_dir)


contrast_strings = differential_expression.make_contrast_strings(comparison_indices, groups_array)

