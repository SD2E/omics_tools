import pandas as pd
import json
import re
import os
import subprocess
import time
from omics_tools import differential_expression, utils, comparison_generator
from collections import Counter

def qc_update(df, factors_to_keep, bool_factors, int_factors):
    df = df[ (df['QC_gcorr_BOOL']==True) & (df['QC_nmap_BOOL']==True) ]

    patterns_to_filter = ["qc_", "_unit", "_input_state"]
    columns_to_filter = ["replicate", "sample_id", "temperature", "timepoint"]
    for col in df.columns:
        if col not in factors_to_keep and (any(p in col.lower() for p in patterns_to_filter) or col.lower() in columns_to_filter):
            df.drop(col, axis=1, inplace=True)
    
    strain_column = ""
    for col in bool_factors:
        df[col] = df[col].astype(bool)
        
    for col in int_factors:
        df[col] = pd.to_numeric(df[col]).astype('int64')

    return df

def create_additive_design(df,int_cols=['IPTG','arabinose']):
    #genes = set(df.index)
    #genes_index = pd.DataFrame(list(range(len(genes))),columns=['Num_Index'],index=genes)
    #df = df.join(genes_index)
    for col in df.columns:
        if col in int_cols:
            df[col]=df[col].astype(int)
    df = pd.get_dummies(df,columns=['Timepoint'])
    #df_test = pd.get_dummies(df,columns=['Num_Index'])
    #return df_test
    return df
    
def main(counts_df_path, result_dir):
    if "29422" in counts_df_path:
        nand20 = False
    else:
        nand20 = True
    counts_df = pd.read_csv(counts_df_path, sep=',', low_memory=False)
    counts_df.rename({counts_df.columns[0]:'sample_id'},inplace=True,axis=1)

    print(counts_df.shape)
    base_factor = ['Strain']

    if nand20:
        # for 23299
        bool_factors=['IPTG','Arabinose']
        local_test = True
        if local_test:
            DE_tests = [
                ['MG1655','MG1655_LPV3'],
                ['MG1655','MG1655_LPV3_AraC_Sensor'],
                ['MG1655','MG1655_LPV3_AraC_Sensor_pBADmin_YFP'],
                ['MG1655','MG1655_LPV3_LacI_AraC_Sensors'],
                ['MG1655','MG1655_LPV3_LacI_AraC_Sensors_pBADmin_PhlF'],
                ['MG1655','MG1655_LPV3_LacI_AraC_Sensors_pBADmin_PhlF_pPhlF_YFP']
            ]
        else:
            DE_tests = [
                ['MG1655','MG1655_LPV3'],
                ['MG1655','MG1655_LPV3_AraC_Sensor'],
                ['MG1655','MG1655_LPV3_AraC_Sensor_pBADmin_YFP'],
                ['MG1655','MG1655_LPV3_LacI_AraC_Sensors'],
                ['MG1655','MG1655_LPV3_LacI_AraC_Sensors_pBADmin_PhlF'],
                ['MG1655','MG1655_LPV3_LacI_AraC_Sensors_pBADmin_PhlF_pPhlF_YFP'],
                ['MG1655','MG1655_LPV3_LacI_AraC_Sensors_pBADmin_PhlF_pTac_AmeR'],
                ['MG1655','MG1655_LPV3_LacI_AraC_Sensors_pBADmin_PhlF_pTac_AmeR_pPhlF_pAmeR_YFP'],
                ['MG1655','MG1655_LPV3_LacI_AraC_Sensors_pBADmin_PhlF_pTac_BM3R1'],
                ['MG1655','MG1655_LPV3_LacI_AraC_Sensors_pBADmin_PhlF_pTac_BM3R1_pPhlF_pBM3R1_YFP'],
                ['MG1655','MG1655_LPV3_LacI_AraC_Sensors_pBADmin_PsrA'],
                ['MG1655','MG1655_LPV3_LacI_AraC_Sensors_pBADmin_PsrA_pPsrA_YFP'],
                ['MG1655','MG1655_LPV3_LacI_AraC_Sensors_pBADmin_PsrA_pTac_AmeR'],
                ['MG1655','MG1655_LPV3_LacI_AraC_Sensors_pBADmin_PsrA_pTac_AmeR_pPsrA_pAmeR_YFP'],
                ['MG1655','MG1655_LPV3_LacI_AraC_Sensors_pBADmin_YFP'],
                ['MG1655','MG1655_LPV3_LacI_AraC_Sensors_pTac_AmeR'],
                ['MG1655','MG1655_LPV3_LacI_AraC_Sensors_pTac_AmeR_pAmeR_YFP'],
                ['MG1655','MG1655_LPV3_LacI_AraC_Sensors_pTac_BM3R1'],
                ['MG1655','MG1655_LPV3_LacI_AraC_Sensors_pTac_BM3R1_pBM3R1_YFP'],
                ['MG1655','MG1655_LPV3_LacI_AraC_Sensors_pTac_YFP'],
                ['MG1655','MG1655_LPV3_LacI_Sensor'],
                ['MG1655','MG1655_LPV3_LacI_Sensor_pTac_AmeR'],
                ['MG1655','MG1655_LPV3_LacI_Sensor_pTac_AmeR_pAmeR_YFP'],
                ['MG1655','MG1655_LPV3_LacI_Sensor_pTac_BM3R1'],
                ['MG1655','MG1655_LPV3_LacI_Sensor_pTac_BM3R1_pBM3R1_YFP'],
                ['MG1655','MG1655_LPV3_LacI_Sensor_pTac_PhlF'],
                ['MG1655','MG1655_LPV3_LacI_Sensor_pTac_PhlF_pPhlF_YFP'],
                ['MG1655','MG1655_LPV3_LacI_Sensor_pTac_PsrA'],
                ['MG1655','MG1655_LPV3_LacI_Sensor_pTac_PsrA_pPsrA_YFP'],
                ['MG1655','MG1655_LPV3_LacI_Sensor_pTac_YFP']
            ]
    else:
        # for 29422
        bool_factors = ['IPTG', 'Cuminic_acid', 'Vanillic_acid', 'Xylose']
        DE_tests = [['Bacillus subtilis 168 Marburg', 'Bacillus subtilis 168 Marburg']]
            
    int_factors = ['Timepoint']
    sub_factors = bool_factors + int_factors
    factors_to_keep = base_factor + sub_factors
    print("factors_to_keep: {}".format(factors_to_keep))

    control_factors = {}
    for i in int_factors:
        if nand20:
            control_factors[i] = 5
        else:
            control_factors[i] = 0
    for bf in bool_factors:
        control_factors[bf] = False
    print("control_factors: {}".format(control_factors))

    counts_df_qcd = qc_update(counts_df, factors_to_keep, bool_factors, int_factors)
    print(counts_df_qcd.shape)
    print(counts_df_qcd.head(5))
    counts_df_qcd.reset_index(inplace=True,drop=True)
    print("Unique strains: {}".format(len(counts_df_qcd['Strain'].unique())))

    run_dir = os.path.join(os.getcwd(), result_dir)
    if not os.path.exists(run_dir):
        os.makedirs(run_dir)
    
    #groups_array = utils.group_by_factors(counts_df_qcd, factors_to_keep)
    
    # This generates aggregate_dict.json
    comparison_indices = comparison_generator.generate_comparisons(counts_df_qcd,
                                                                   base_comparisons = DE_tests,
                                                                   base_factor = base_factor, 
                                                                   sub_factors = sub_factors,
                                                                   freedom = len(sub_factors),
                                                                   aggregation_flag = True,
                                                                   run_dir = run_dir,
                                                                   control_factor_in = control_factors)
    
    #contrast_strings = differential_expression.make_contrast_strings(comparison_indices, groups_array)
    
    #print("Counter(groups_array): {}".format(Counter(groups_array)))
    print("len(comparison_indices): {}".format(len(comparison_indices)))
    
    # The code below shows two ways of launching R, only the scripts approach has been proven
    gen_r_scripts = True
    if gen_r_scripts:
        differential_expression.make_hpc_de_files(dataframe = counts_df_qcd,
                                                  base_comparisons = DE_tests,
                                                  sub_factors = sub_factors,
                                                  run_dir = run_dir,
                                                  filter_unused_base_factors = True,
                                                  freedom = len(sub_factors),
                                                  export_tagwise_noise = False,
                                                  base_factor = base_factor,
                                                  control_factor_in = control_factors)
        os.chdir(result_dir)
        print("new working dir: {}".format(os.getcwd()))
        subprocess.call(['chmod', 'u+x', './dge_local.sh'])
        ret_value = subprocess.call(['./dge_local.sh'])
    else:
        r_cmds = differential_expression.make_DE_cmds(dataframe = counts_df_qcd,
                                                      base_comparisons = DE_tests,
                                                      base_factor = base_factor, 
                                                      sub_factors = sub_factors,
                                                      control_factor_in = control_factors)
        print('Created {0} differential tests'.format(len(r_cmds)))
        print("r_cmds: {}".format(r_cmds))

    print("ret_value: {}".format(ret_value))
    if ret_value == 0:
        while True:
            result_files = [f for f in os.listdir('./results') if re.match(r'^(?!dge).*$', f)]
            if len(result_files) < len(comparison_indices):
                time.sleep(10)
            else:
                break

        df_diff_exp = utils.aggregate_dataframes(run_dir, sub_factors, suffix_col='Strain_2')
        dfs=[]
        for test in DE_tests:
            strain_2 = test[1]
            cols_of_interest = ['FDR_'+strain_2,'nlogFDR_'+strain_2,'logFC_'+strain_2]+sub_factors
            df_sub = df_diff_exp[cols_of_interest]
            df_sub.rename({'FDR_'+strain_2:'FDR','nlogFDR_'+strain_2:'nlogFDR','logFC_'+strain_2:'logFC'},axis=1,inplace=True)
            df_sub['strain']=strain_2
            dfs.append(df_sub)
        df_diff_exp_all = pd.concat(dfs)
        df_diff_additive_design = create_additive_design(df_diff_exp_all, int_cols=bool_factors)
        print("df_diff_additive_design")
        print(df_diff_additive_design.head(5))
        df_diff_additive_design.to_csv(os.path.join(run_dir,'results/additive_design_df.csv'))
        
if __name__ == '__main__':
    #main("./experiment.ginkgo.29422_ReadCountMatrix_preCAD_transposed.csv", '../exp_ref_additive_design')
    main("./experiment.ginkgo.23299_ReadCountMatrix_preCAD_transposed.csv", '../exp_ref_additive_design')