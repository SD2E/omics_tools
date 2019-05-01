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