import argparse
import pandas as pd
import seaborn as sns
import json
import re
import os
import copy
import subprocess
import time
from omics_tools import differential_expression, utils, comparison_generator
from collections import Counter, defaultdict

from multiprocessing import Pool, get_context
import parallel_comparison

import matplotlib.pyplot as plt

# requests prefers simplejson
try:
    import simplejson as json
    from simplejson.errors import JSONDecodeError
except ImportError:
    import json
    from json.decoder import JSONDecodeError

#%matplotlib inline
#import seaborn as sns
#sns.set()

def qc_update(df, factors_to_keep, bool_factors, int_factors):
    df = df[ (df['qc_gcorr_bool']==True) & (df['qc_nmap_bool']==True) ]

    patterns_to_filter = ["qc_", "_unit", "_input_state"]
    columns_to_filter = ["replicate", "sample_id", "temperature", "timepoint","dextrose"]
    for col in df.columns:
        if col not in factors_to_keep and (any(p in col.lower() for p in patterns_to_filter) or col.lower() in columns_to_filter):
            df.drop(col, axis=1, inplace=True)
    
    strain_column = ""
    for col in bool_factors:
        df[col] = df[col].astype(bool)
        
    for col in int_factors:
        df[col] = pd.to_numeric(df[col]).astype('int64')

    return df

def create_additive_design(df,int_cols=['iptg','arabinose']):
    #genes = set(df.index)
    #genes_index = pd.DataFrame(list(range(len(genes))),columns=['Num_Index'],index=genes)
    #df = df.join(genes_index)
    for col in df.columns:
        if col in int_cols:
            df[col]=df[col].astype(int)
    #df = pd.get_dummies(df,columns=['Timepoint'])
    #df_test = pd.get_dummies(df,columns=['Num_Index'])
    #return df_test
    return df

def comparison_heatmap(cfm_input_df, exp_condition_cols, target_column, replicates=True, figure_name='comparison_heatmap'):
    additional_conditions = []
    if replicates == True:
        additional_conditions.append('replicate')
        
    print("exp_condition_cols: {}".format(exp_condition_cols))
    
    start = time.perf_counter()
    prior_comparisons = set()
    execution_space = list()
    data = defaultdict(list)
    condition_groups = cfm_input_df.groupby(exp_condition_cols + additional_conditions)
    print("cfm_input_df")
    print(cfm_input_df.head(5))
    cols_of_interest0 = cfm_input_df[exp_condition_cols + additional_conditions].drop_duplicates()
    sort_keys = exp_condition_cols + additional_conditions
    cols_of_interest1 = cols_of_interest0.sort_values(sort_keys)
    cols_of_interest = cols_of_interest1.values

    for a,condition1 in enumerate(cols_of_interest):
        for b,condition2 in enumerate(cols_of_interest):
#             current_comparison = frozenset(Counter((a,b)))
            if (a,b) not in prior_comparisons:
                execution_space.append((condition1, condition2, condition_groups, target_column))
                prior_comparisons.add((a,b))
                prior_comparisons.add((b,a))

    num_processors=os.cpu_count()
    pool = Pool(processes = num_processors)
    output = pool.map(parallel_comparison.perform_matrix_calculation, execution_space)
    pool.close()

    for comparison_i,item in enumerate(execution_space):
        condition1,condition2,_pass1,_pass2 = item
        data['condition1'].append(", ".join(map(str, condition1)))
        data['condition2'].append(", ".join(map(str, condition2)))
        for i,variable in enumerate(condition1):
            data['{}_1'.format((exp_condition_cols + additional_conditions)[i])].append(condition1[i])
        for i,variable in enumerate(condition2):
            data['{}_2'.format((exp_condition_cols + additional_conditions)[i])].append(condition2[i])
        data['comparison'].append(output[comparison_i])
        data['condition1'].append(", ".join(map(str, condition2)))
        data['condition2'].append(", ".join(map(str, condition1)))
        for i,variable in enumerate(condition1):
            data['{}_1'.format((exp_condition_cols + additional_conditions)[i])].append(condition1[i])
        for i,variable in enumerate(condition2):
            data['{}_2'.format((exp_condition_cols + additional_conditions)[i])].append(condition2[i])
        data['comparison'].append(output[comparison_i])

    df = pd.DataFrame.from_dict(data)
    
    df = df.drop_duplicates()

    pivoted_df = df.pivot('condition1','condition2','comparison')

    from matplotlib.colors import LinearSegmentedColormap
    plt.figure(figsize=(60, 60))
    sns.heatmap(pivoted_df,linecolor='black',cmap=LinearSegmentedColormap.from_list('try1', ['red', 'white'], N=10))
    b,t = plt.ylim()
    b += 0.5
    t -= 0.5
    plt.ylim(b,t)
    plt.tight_layout()
    # plt.show()
    plt.savefig(figure_name + '.png')

    stop = time.perf_counter()
    print(stop-start)
    return df

def load_config(config_file):
    try:
        with open(config_file) as json_data:
            config_json = json.load(json_data)
    except Exception as e:
        print('Failed to load config_file', e)

    int_factors = config_json["int_factors"] if "int_factors" in config_json else []
    bool_factors = config_json["bool_factors"] if "bool_factors" in config_json else []
    float_factors = config_json["float_factors"] if "float_factors" in config_json else []
    other_factors = config_json["other_factors"] if "other_factors" in config_json else []
    parts = config_json["parts"] if "parts" in config_json else []
    hrm_experimental_condition_cols = int_factors + bool_factors + float_factors + other_factors + parts
    control_factors = config_json["control_factors"] if "control_factors" in config_json else {}
    cf_value = config_json["cf_value"]
    DE_tests = config_json["DE_tests"]
    add_one = config_json["add_one"]
    fdr_max = config_json["fdr_max"]
    log_fc_min = config_json["log_fc_min"]
    batch_delay = config_json["batch_delay"]
    output_name = config_json["output_name"]
    comparison_target_column = config_json["comparison_target_column"]
    if 'base_factor' in config_json:
        base_factor = [config_json['base_factor']]
    else:
        base_factor = None
    return int_factors, bool_factors, float_factors, control_factors, cf_value, hrm_experimental_condition_cols, DE_tests, add_one, fdr_max, log_fc_min, batch_delay, output_name, comparison_target_column,base_factor

def main(counts_df_path, config_file, result_dir):

    counts_df = pd.read_csv(counts_df_path, sep=',', low_memory=False)
    counts_df.rename({counts_df.columns[0]:'sample_id'},inplace=True,axis=1)
    counts_df.rename({col: str.lower(col) for col in counts_df.columns},inplace=True, axis=1)

    print(counts_df.shape)


    int_factors, bool_factors, float_factors, control_factors, cf_value, hrm_experimental_condition_cols, DE_tests, add_one, fdr_max, log_fc_min, batch_delay, output_name, comparison_target_column,base_factor = load_config(config_file)
    print("hrm_experimental_condition_cols: {}".format(hrm_experimental_condition_cols))
    if not base_factor:
        base_factor = ['strain']
    sub_factors = int_factors + bool_factors + float_factors 
    factors_to_keep = base_factor + sub_factors
    print("factors_to_keep: {}".format(factors_to_keep))

    for i in int_factors:
        if i not in control_factors:
            control_factors[i] = cf_value

    for bf in bool_factors:
        control_factors[bf] = False
    # for ff in float_factors:
    #     control_factors[ff] = 0
        
    print("control_factors: {}".format(control_factors))

    counts_df_qcd = qc_update(counts_df, factors_to_keep, bool_factors, int_factors)
    print(counts_df_qcd.shape)
    print(counts_df_qcd.head(5))
    counts_df_qcd.reset_index(inplace=True,drop=True)
    try:
        print("Unique strains: {}".format(len(counts_df_qcd['strain'].unique())))
    except:
        print('No column named strain')

    run_dir = os.path.join(os.getcwd(), result_dir)
    if not os.path.exists(run_dir):
        os.makedirs(run_dir)
        
    # May need to remove content of the result_dir
    
    #groups_array = utils.group_by_factors(counts_df_qcd, factors_to_keep)
    
    # This generates aggregate_dict.json
    comparison_indices = comparison_generator.generate_comparisons(counts_df_qcd,
                                                                   base_comparisons = DE_tests,
                                                                   base_factor = base_factor, 
                                                                   sub_factors = sub_factors,
                                                                   freedom = len(sub_factors)+1 if add_one else len(sub_factors),
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
#                                                  aggregation_flag = True,
                                                  base_comparisons = DE_tests,
                                                  sub_factors = sub_factors,
                                                  run_dir = run_dir,
                                                  filter_unused_base_factors = True,
                                                  freedom = len(sub_factors)+1 if add_one else len(sub_factors),
                                                  export_tagwise_noise = False,
                                                  base_factor = base_factor,
                                                  control_factor_in = control_factors,
                                                  batch_delay=batch_delay)

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
                #time.sleep(10)
                pass
            else:
                break

        df_diff_exp = utils.aggregate_dataframes(run_dir, sub_factors, suffix_col=base_factor[0]+'_2')
        dfs=[]
        for test in DE_tests:
            strain_2 = test[1]
            cols_of_interest = ['flagedgeRremoved_'+strain_2,'FDR_'+strain_2,'nlogFDR_'+strain_2,'logFC_'+strain_2]+sub_factors
            df_sub = df_diff_exp[cols_of_interest]
            df_sub.rename({'FDR_'+strain_2:'FDR','nlogFDR_'+strain_2:'nlogFDR','logFC_'+strain_2:'logFC'},axis=1,inplace=True)
            df_sub['strain']=strain_2
            dfs.append(df_sub)
        df_diff_exp_all = pd.concat(dfs)
        print("df_diff_exp_all")
        print(df_diff_exp_all.head(5))
        df_diff_additive_design = create_additive_design(df_diff_exp_all, int_cols=bool_factors)
        print("df_diff_additive_design")
        print(df_diff_additive_design.head(5))
        df_diff_additive_design.to_csv(os.path.join(run_dir,'results/additive_design_df.csv'))

        hrm_data = df_diff_additive_design.copy()
        hrm_data["gene"] = hrm_data.index
        print("after adding index")
        print(hrm_data.head(5))

        #hrm_data.to_csv("hrm_data_before.csv")

        hrm_data = hrm_data[(hrm_data['FDR'] < fdr_max) & (hrm_data['logFC'] > log_fc_min)]
        #hrm_data.to_csv("hrm_data_after.csv")
        output = comparison_heatmap(hrm_data,hrm_experimental_condition_cols,comparison_target_column,replicates=False,figure_name=output_name)
        output.to_csv(output_name + "_overlap.csv")
        print("output")
        print(output.head(5))
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_counts_file", help="input file")
    parser.add_argument("--config_file", help="analysis configuration")
    parser.add_argument("--output_dir", help="results folder")

    args = parser.parse_args()
    arg_input_counts_file = args.input_counts_file
    arg_config_file = args.config_file
    arg_output_dir = args.output_dir
    
    main(arg_input_counts_file, arg_config_file, arg_output_dir)
    
    #main("./experiment.ginkgo.29422_ReadCountMatrix_preCAD_transposed.csv", '../exp_ref_additive_design')
    #main("./experiment.ginkgo.23299_ReadCountMatrix_preCAD_transposed.csv", '../exp_ref_additive_design')