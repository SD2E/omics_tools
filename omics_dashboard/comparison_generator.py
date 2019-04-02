from itertools import product
from functools import reduce
import operator


def generate_comparisons(df, base_comparisons=None, base_factor=['strain'], sub_factors=None, freedom=1):
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
    for base1, base2 in base_comparisons:
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
        c1_i = df[reduce(operator.and_, ((df[x] == c1[i]) for i, x in enumerate(factors)))].index
        c2_i = df[reduce(operator.and_, ((df[x] == c2[i]) for i, x in enumerate(factors)))].index
        if comp[::-1] in comp_indices:
            continue
        if c1_i.empty or c2_i.empty:
            continue
        # This next test should never happen, we filter it out when constructing the initial DF
        if len(c1_i) == 1 or len(c2_i) == 1:
            continue
        comp_indices[tuple(map(tuple, comp))] = [c1_i, c2_i]
    return comp_indices

    """
    all_comparison_sets = []
    contrast_indices = {}
    for group in groups:
        for contrast in base_comparisons:
            if group[0] in contrast:
                levels = len(factors)
                while levels >= len(factors)-degree:  # 4
                    entry = (group[0], set(factors[:levels]))

                    if levels not in contrast_indices:
                        contrast_indices[levels] = {}

                    group_subset = df.groupby(list(entry[1])).groups
                    for group_in_subset in group_subset:
                        if entry[0] in group_in_subset:
                            if isinstance(group_in_subset, str):
                                comparison = (group_in_subset,)
                                if not set(comparison) in all_comparison_sets:
                                    all_comparison_sets.append(set(comparison))
                            else:
                                comparison = group_in_subset
                                if not set(group_in_subset) in all_comparison_sets:
                                    all_comparison_sets.append(set(group_in_subset))

                            if comparison not in all_comparison_sets:
                                if isinstance(group_in_subset, tuple):
                                    contrast_indices[levels][group_in_subset] = group_subset[group_in_subset]
                                else:
                                    contrast_indices[levels][(group_in_subset,)] = group_subset[group_in_subset]

                    levels -= 1
                    if len(entry[1]) == 1:
                        levels -= 1

                levels = len(factors)
                while levels >= len(factors)-degree:  # 2

                    entry = (group[0], tuple(r_factors[:levels]))
                    group_subset = df.groupby(list(entry[1])).groups
                    for group_in_subset in group_subset:
                        if entry[0] in group_in_subset:

                            if isinstance(group_in_subset, str):
                                comparison = (group_in_subset,)
                                if not set(comparison) in all_comparison_sets:
                                    all_comparison_sets.append(set(comparison))
                                else:
                                    continue
                            else:
                                comparison = group_in_subset
                                if not set(comparison) in all_comparison_sets:
                                    all_comparison_sets.append(set(group_in_subset))
                                else:
                                    continue

                            if comparison not in all_comparison_sets:
                                if isinstance(group_in_subset, tuple):
                                    contrast_indices[levels][group_in_subset] = group_subset[group_in_subset]
                                else:
                                    contrast_indices[levels][(group_in_subset,)] = group_subset[group_in_subset]
                    levels -= 1
                    if len(entry[1]) == 1:
                        levels -= 1

    comparison_indices = {}
    for level in range(len(factors)-degree, len(factors) + 1):
        level_comparisons = []
        for m in contrast_indices[level].keys():
            for n in contrast_indices[level].keys():
                for group in base_comparisons:
                    group1, group2 = group[0], group[1]
                    if any([True for x in m if x == group1]) and any([True for x in n if x == group2]):
                        free = [1 for x in m for y in n if x == y]
                        if sum(free) == level - freedom:
                            cats_m = [factor_to_categories[x.strip()] for x in m]
                            cats_n = [factor_to_categories[x.strip()] for x in n]
                            if cats_m == cats_n:
                                if [m, n][::-1] not in level_comparisons:
                                    level_comparisons.append([m, n])

        for lev in level_comparisons:
            comparison_indices[tuple(lev)] = []
            collect = []
            for comparison in lev:
                for ind in contrast_indices[level]:
                    if ind == comparison:
                        collect.append(contrast_indices[level][ind])
            comparison_indices[tuple(lev)] = collect
    return comp_indices
    """
