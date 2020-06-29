from itertools import product
from functools import reduce
import operator
import os
import json

def generate_comparisons(df, base_comparisons=None, base_factor=['strain'], sub_factors=None, freedom=1,aggregation_flag=False,run_dir=None,control_factor_in=None):
    """
    :param df: pandas.DataFrame with all the experimental metadata and data.
    :param base_comparisons: List of list with two strings declaring base level comparisons.
                             [[Strain1, Strain2], [Strain2, Strain9]
    :param base_factor: List containing the column name containing the identifiers for
                        the base comparisons.
    :param sub_factors: List of factors to use as subdividers of the base comparisons.
                        [time, temp] ~ Strain1_time1_temp1 vs Strain2_time1_temp1
                                     ~ Strain1_time1_temp1 vs Strain2_time1_temp2
                                     ~ Strain1_time1_temp2 vs Strain2_time1_temp1
                                     ~ Strain1_time2_temp1 vs Strain2_time1_temp1
                                       etc.
  :param control_factor_in: Dictionary of the format {subfactor:subfactor_value} to be used as the single control
                        {'timepoint':'time1','temperature':temp1}
                        note that all keys in the dictionary MUST be in the subfactor list.
    :return comparisons_indices:
    """
    comparisons = set()
    sub_factors = sorted(sub_factors)

    control_factor = {}
    for item in sub_factors:
        control_factor[item]=control_factor_in[item]
    factors = base_factor + sub_factors

    if control_factor:
        for key in control_factor:
            if key not in sub_factors:
                raise RuntimeError('{0} not present in factors'.format(key))
        if len(sub_factors)!=len(control_factor.keys()):
            raise RuntimeError('Control factor is not assigning a value to each subfactor')

    agg_dict = {}
    for base1, base2 in base_comparisons:
        if control_factor:
            factor_enums = [set([(control_factor[x])]) for x in control_factor]

        else:
            factor_enums = [set(df[df[base_factor[0]] == base1][x]) for x in sub_factors]

        comparisons_1 = list(product(*factor_enums))
        factor_enums = [set(df[df[base_factor[0]] == base2][x]) for x in sub_factors]
        comparisons_2 = list(product(*factor_enums))


        for i in range(len(factors) - freedom, len(factors)):
            valids = []
            for comp1 in comparisons_1:
                for comp2 in comparisons_2:
                    same = 0
                    for j in range(len(sub_factors)):
                        if comp1[j] == comp2[j]:
                            same += 1
                    if base1 == base2:
                        if same == (i - 1):
                            valids.append([[base1] + list(comp1), [base2] + list(comp2)])
                    else:
                        if same == i:
                            valids.append([[base1] + list(comp1), [base2] + list(comp2)])
            if valids:
                for v in valids:
                    comparisons.add(tuple(map(tuple, v)))

    comp_indices = {}
    for comp in comparisons:
        c1, c2 = comp[0], comp[1]
        c1_str = '_'.join(map(lambda x: str(x), c1))
        c2_str = '_'.join(map(lambda x: str(x), c2))
        fname = c1_str+'-vs-'+c2_str+'.txt'

        #create aggregated dataframe dictionary that will be written out to file
        #dictionary is of the form {filename:dict(metadata)}
        agg_dict[fname]={}

        for kk in range(len(c1)):
            agg_dict[fname][factors[kk]+'_1'] = c1[kk]
        for kk in range(len(c2)):
            agg_dict[fname][factors[kk]+'_2'] = c2[kk]

        c1_i = df[reduce(operator.and_, ((df[x] == c1[i]) for i, x in enumerate(factors)))].index
        c2_i = df[reduce(operator.and_, ((df[x] == c2[i]) for i, x in enumerate(factors)))].index
        if comp[::-1] in comp_indices:
            continue
        if c1_i.empty or c2_i.empty:
            print('Removing test', c1_str, 'vs', c2_str, 'because one of them is absent in the data')
            continue
        # This next test should never happen, we filter it out when constructing the initial DF
        if len(c1_i) == 1 or len(c2_i) == 1:
            print('Removing test',c1_str,'vs',c2_str,'because there is only a single replicate')
            continue

        comp_indices[tuple(map(tuple, comp))] = [c1_i, c2_i]

    if aggregation_flag:
        if run_dir is None:
            run_dir = os.getcwd()
        if not os.path.exists(run_dir):
            os.mkdir(run_dir)
        if not os.path.exists(os.path.join(run_dir,'tmp')):
            os.mkdir(os.path.join(run_dir,'tmp'))
            with open(os.path.join(run_dir,'tmp','aggregate_dict.json'), 'w+') as fp:
                json.dump(agg_dict, fp)

    return comp_indices


def generate_comparisons_with_single_control(df, base_comparisons=None, base_factor=['strain'], sub_factors=None, freedom=1,aggregation_flag=False,run_dir=None,control_factor=None):
    """
    :param df: pandas.DataFrame with all the experimental metadata and data.
    :param base_comparisons: List of list with two strings declaring base level comparisons.

                             [[Strain1, Strain2], [Strain2, Strain9]
    :param base_factor: List containing the column name containing the identifiers for
                        the base comparisons.
    :param sub_factors: List of factors to use as subdividers of the base comparisons.
                        [time, temp] ~ Strain1_time1_temp1 vs Strain2_time1_temp1
                                     ~ Strain1_time1_temp1 vs Strain2_time1_temp2
                                     ~ Strain1_time1_temp2 vs Strain2_time1_temp1
                                     ~ Strain1_time2_temp1 vs Strain2_time1_temp1
                                       etc.

    :return comparisons_indices:
    """
    comparisons = set()
    sub_factors = sorted(sub_factors)
    factors = base_factor + sub_factors

    agg_dict = {}
    for base1, base2 in base_comparisons:
        factor_enums = [set(control_factor[x]) for x in control_factor]
        comparisons_1 = list(product(*factor_enums))
        factor_enums = [set(df[df[base_factor[0]] == base2][x]) for x in sub_factors]
        comparisons_2 = list(product(*factor_enums))


        for i in range(len(factors) - freedom, len(factors)):
            valids = []
            for comp1 in comparisons_1:
                for comp2 in comparisons_2:
                    same = 0
                    for j in range(len(sub_factors)):
                        if comp1[j] == comp2[j]:
                            same += 1
                    if base1 == base2:
                        if same == (i - 1):
                            valids.append([[base1] + list(comp1), [base2] + list(comp2)])
                    else:
                        if same == i:
                            valids.append([[base1] + list(comp1), [base2] + list(comp2)])
            if valids:
                for v in valids:
                    comparisons.add(tuple(map(tuple, v)))

    comp_indices = {}
    for comp in comparisons:
        c1, c2 = comp[0], comp[1]
        c1_str = '_'.join(map(lambda x: str(x), c1))
        c2_str = '_'.join(map(lambda x: str(x), c2))
        fname = c1_str+'-vs-'+c2_str+'.txt'

        #create aggregated dataframe dictionary that will be written out to file
        #dictionary is of the form {filename:dict(metadata)}
        agg_dict[fname]={}

        for kk in range(len(c1)):
            agg_dict[fname][factors[kk]+'_1'] = c1[kk]
        for kk in range(len(c2)):
            agg_dict[fname][factors[kk]+'_2'] = c2[kk]

        c1_i = df[reduce(operator.and_, ((df[x] == c1[i]) for i, x in enumerate(factors)))].index
        c2_i = df[reduce(operator.and_, ((df[x] == c2[i]) for i, x in enumerate(factors)))].index
        if comp[::-1] in comp_indices:
            continue
        if c1_i.empty or c2_i.empty:
            print('Removing test', c1_str, 'vs', c2_str, 'because one of them is absent in the data')
            continue
        # This next test should never happen, we filter it out when constructing the initial DF
        if len(c1_i) == 1 or len(c2_i) == 1:
            print('Removing test',c1_str,'vs',c2_str,'because there is only a single replicate')
            continue

        comp_indices[tuple(map(tuple, comp))] = [c1_i, c2_i]

    if aggregation_flag:
        if run_dir is None:
            run_dir = os.getcwd()
        if not os.path.exists(run_dir):
            os.mkdir(run_dir)
        if not os.path.exists(os.path.join(run_dir,'tmp')):
            os.mkdir(os.path.join(run_dir,'tmp'))
            with open(os.path.join(run_dir,'tmp','aggregate_dict.json'), 'w+') as fp:
                json.dump(agg_dict, fp)

    return comp_indices

