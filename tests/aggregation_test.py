from omics_tools.utils import aggregate_dataframes


run_dir = '/Users/meslami/Desktop/DARPA/SD2/NovelChassis/NAND_20_additive_design'
sub_factors = ['timepoint', 'IPTG', 'arabinose']
aggregate_dataframes(run_dir,sub_factors,suffix_col='strain_2')