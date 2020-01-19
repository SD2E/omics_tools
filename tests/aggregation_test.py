from omics_tools.utils import aggregate_dataframes


run_dir = '/Users/meslami/Desktop/DARPA/SD2/NovelChassis/Bacillus/Inducer_Prediction'
sub_factors = ['timepoint']
aggregate_dataframes(run_dir,sub_factors,suffix_col='strain_inducer_2')