from scipy import stats


def perform_matrix_calculation(condition_grouping,target_col='gene'):
	condition1,condition2,grouped_df = condition_grouping
	comparison = 0
	if target_col == 'BL1-A_MEFL':
		comparison = compute_EMD(grouped_df, condition1, condition2, target_col)
	else:
		comparison = compute_set_overlap(grouped_df, condition1, condition2, target_col)
	return comparison

def retrieve_values(grouped_df,desired_condition,target_col='gene'):
	desired_condition = tuple(desired_condition)
	temp_df = grouped_df.get_group(desired_condition)
	return temp_df[target_col].values

def compute_EMD(grouped_df, desired_condition1, desired_condition2, target_col='BL1-A_MEFL'):
	values1 = retrieve_values(grouped_df, desired_condition1, target_col)
	values2 = retrieve_values(grouped_df, desired_condition2, target_col)
	EMD  = stats.wasserstein_distance(values2, values1)
	return EMD

def compute_set_overlap(grouped_df, desired_condition1, desired_condition2, target_col='gene'):
	values1 = retrieve_values(grouped_df, desired_condition1, target_col)
	values2 = retrieve_values(grouped_df, desired_condition2, target_col)
	
	s1 = set(values2)
	s2 = set(values1)
	
	return float(len(s1.intersection(s2))) / float(len(s1.union(s2)))