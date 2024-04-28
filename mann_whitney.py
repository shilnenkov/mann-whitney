import pandas as pd
import os
import re
from scipy.stats import mannwhitneyu, chi2_contingency, fisher_exact

SEPARATOR = ';'
CSV_PATH = 'data/oncolab_freq/list_2.csv'


FLD_EGFR_TYPE = 'EGFR type'
FLD_AGE = 'Возр'
FLD_SEX = 'Пол'
FLD_SMOKING = 'Статус курения'


def read_csv(filepath: str) -> pd.DataFrame:
    data = pd.read_csv(filepath, sep=SEPARATOR)
    return data[data[FLD_AGE] > 1]


def log_to_file(data: pd.DataFrame, sheet_name: str) -> None:
    OUT_DIR = 'логи'

    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)
    with pd.ExcelWriter(f'{OUT_DIR}/{sheet_name}.xlsx', engine='openpyxl', mode='w') as writer:
        data.to_excel(writer, index=False)


def mutation_destrib_by_ages(data: pd.DataFrame, min_age: int, max_age: int, title: str) -> None:
    is_greater_than_min_age = data[FLD_AGE] > min_age
    is_less_or_eq_max_age = data[FLD_AGE] <= max_age
    age_condition = is_greater_than_min_age & is_less_or_eq_max_age

    all = data[FLD_AGE].size
    count = data[age_condition][FLD_AGE].size
    percent = count / all * 100.0

    print(f'{title}: {min_age}-{max_age}: {count}/{all} ({percent:.2f}%)')


def compare_with_and_without_mutations(data: pd.DataFrame) -> None:
    has_no_mutations_condition = data[FLD_EGFR_TYPE].str.contains('WT')
    has_mutations_condition = ~(has_no_mutations_condition)

    with_mutations = data[has_mutations_condition][FLD_AGE]
    without_mutations = data[has_no_mutations_condition][FLD_AGE]

    statistics, p_value = mannwhitneyu(with_mutations, without_mutations)
    print(f'С мутациями ({with_mutations.size}) и '
          f'без ({without_mutations.size}) - '
          f'Mann-Whitney U statistic: {statistics}, P-value: {p_value}')

    log_to_file(data[has_mutations_condition], 'с_мутациями')
    log_to_file(data[has_no_mutations_condition], 'буз_мутаций')


def compare_freq_and_rare_mutations(data: pd.DataFrame) -> None:
    has_mutations_condition = ~(data[FLD_EGFR_TYPE].str.contains('WT'))

    is_ex19del = data[FLD_EGFR_TYPE].str.contains('ex19del')
    is_l858r = data[FLD_EGFR_TYPE].str.contains('L858R')

    is_s768i = data[FLD_EGFR_TYPE].str.contains('S768I')
    is_l703v = data[FLD_EGFR_TYPE].str.contains('L703V')
    is_g779c = data[FLD_EGFR_TYPE].str.contains('G779C')
    is_d761y = data[FLD_EGFR_TYPE].str.contains('D761Y')

    freq_mutations_condition = (is_ex19del | is_l858r) & \
                               ~(is_s768i | is_l703v | is_g779c | is_d761y)

    is_double_mut = data[FLD_EGFR_TYPE].str.contains(re.escape('+'))

    has_freq_mutations = freq_mutations_condition
    has_rare_mutations = has_mutations_condition & ~(freq_mutations_condition)
    has_rare_double_mutations = has_rare_mutations & is_double_mut

    freq_mutations = data[has_freq_mutations][FLD_AGE]
    rare_mutations = data[has_rare_mutations][FLD_AGE]
    rare_double_mutations = data[has_rare_double_mutations][FLD_AGE]

    statistics, p_value = mannwhitneyu(freq_mutations, rare_mutations)
    print(f'Частые ({freq_mutations.size}) и '
          f'редкие ({rare_mutations.size}) - '
          f'Mann-Whitney U statistic: {statistics}, P-value: {p_value}')

    statistics, p_value = mannwhitneyu(freq_mutations, rare_double_mutations)
    print(f'Частые ({freq_mutations.size}) и '
          f'редкие двойные ({rare_double_mutations.size}) - '
          f'Mann-Whitney U statistic: {statistics}, P-value: {p_value}')

    log_to_file(data[has_freq_mutations], 'частые_мутации')
    log_to_file(data[has_rare_mutations], 'редкие_мутации')
    log_to_file(data[has_rare_double_mutations], 'двойные_мутации')


def calc_fisher_by_sex(data: pd.DataFrame) -> None:
    has_mutation = data[FLD_EGFR_TYPE].str.contains('WT')
    has_no_mutation = ~has_mutation

    men = data[FLD_SEX].str.contains('м')
    women = ~men

    men_have_mut = data[men & has_mutation].shape[0]
    men_have_no_mut = data[men & has_no_mutation].shape[0]
    women_have_mut = data[women & has_mutation].shape[0]
    women_have_no_mut = data[women & has_no_mutation].shape[0]

    print('\t\tБез мутаций\tС мутацией')
    print(f'Мужчины\t{men_have_no_mut}\t\t\t{men_have_mut}')
    print(f'Женщины\t{women_have_no_mut}\t\t\t{women_have_mut}')

    table = [
        [men_have_mut, men_have_no_mut],
        [women_have_mut, women_have_no_mut],
    ]

    odds_ratio, p_value = fisher_exact(table)
    print(f'Фишер: Odds ratio: {odds_ratio}, p_value: {p_value}')

    chi2_statistic, p_value, dof, expected = chi2_contingency(table)
    print(f'chi^2 stat: {chi2_statistic}, p_value: {p_value}')


def calc_fisher_by_smoking(data: pd.DataFrame) -> None:
    known_smokers = data[data[FLD_SMOKING].str.len() > 0]

    not_smokers = known_smokers[FLD_SMOKING].str.contains('не ')
    smokers = ~not_smokers

    has_no_mutation = known_smokers[FLD_EGFR_TYPE].str.contains('WT')
    has_mutation = ~has_no_mutation

    smokers_have_mut = known_smokers[smokers & has_mutation].shape[0]
    smokers_have_no_mut = known_smokers[smokers & has_no_mutation].shape[0]
    not_smokers_have_mut = known_smokers[not_smokers & has_mutation].shape[0]
    not_smokers_have_no_mut = known_smokers[not_smokers & has_no_mutation].shape[0]

    print('\t\t\tБез мутаций\tС мутацией')
    print(f'Курят\t\t{smokers_have_mut}\t\t\t{smokers_have_no_mut}')
    print(f'Не курят\t{not_smokers_have_mut}\t\t\t{not_smokers_have_no_mut}')

    table = [
        [smokers_have_mut, smokers_have_no_mut],
        [not_smokers_have_mut, not_smokers_have_no_mut],
    ]

    odds_ratio, p_value = fisher_exact(table)
    print(f'Odds ratio: {odds_ratio}, p_value: {p_value}')

    chi2_statistic, p_value, dof, expected = chi2_contingency(table)
    print(f'chi^2 stat: {chi2_statistic}, p_value: {p_value}')


def calc_fisher_freq_rare_by_sex(data: pd.DataFrame) -> None:
    has_mutations_condition = ~(data[FLD_EGFR_TYPE].str.contains('WT'))

    is_ex19del = data[FLD_EGFR_TYPE].str.contains('ex19del')
    is_l858r = data[FLD_EGFR_TYPE].str.contains('L858R')

    is_s768i = data[FLD_EGFR_TYPE].str.contains('S768I')
    is_l703v = data[FLD_EGFR_TYPE].str.contains('L703V')
    is_g779c = data[FLD_EGFR_TYPE].str.contains('G779C')
    is_d761y = data[FLD_EGFR_TYPE].str.contains('D761Y')

    freq_mutations_condition = (is_ex19del | is_l858r) & \
                               ~(is_s768i | is_l703v | is_g779c | is_d761y)

    is_double_mut = data[FLD_EGFR_TYPE].str.contains(re.escape('+'))

    has_freq_mutations = freq_mutations_condition
    has_rare_mutations = has_mutations_condition & ~(freq_mutations_condition)
    has_rare_double_mutations = has_rare_mutations & is_double_mut

    men = data[FLD_SEX].str.contains('м')
    women = ~men

    men_have_freq = data[men & has_freq_mutations].shape[0]
    men_have_rare = data[men & has_rare_mutations].shape[0]
    men_have_double_rare = data[men & has_rare_double_mutations].shape[0]
    women_have_freq = data[women & has_freq_mutations].shape[0]
    women_have_rare = data[women & has_rare_mutations].shape[0]
    women_have_double_rare = data[women & has_rare_double_mutations].shape[0]

    print('\t\tРедкие мутации\tЧастые мутации')
    print(f'Мужчины\t{men_have_rare}\t\t\t\t{men_have_freq}')
    print(f'Женщины\t{women_have_rare}\t\t\t\t{women_have_freq}')

    # Freq & double rare

    table = [
        [men_have_rare, men_have_freq],
        [women_have_rare, women_have_freq],
    ]

    odds_ratio, p_value = fisher_exact(table)
    print(f'Фишер: Odds ratio: {odds_ratio}, p_value: {p_value}')

    chi2_statistic, p_value, dof, expected = chi2_contingency(table)
    print(f'chi^2 stat: {chi2_statistic}, p_value: {p_value}')

    print('\t\tЧастые мутации\tРедкие двойные мутации')
    print(f'Мужчины\t{men_have_freq}\t\t\t\t{men_have_double_rare}')
    print(f'Женщины\t{women_have_freq}\t\t\t\t{women_have_double_rare}')

    table = [
        [men_have_freq, men_have_double_rare],
        [women_have_freq, women_have_double_rare],
    ]

    odds_ratio, p_value = fisher_exact(table)
    print(f'Фишер: Odds ratio: {odds_ratio}, p_value: {p_value}')

    chi2_statistic, p_value, dof, expected = chi2_contingency(table)
    print(f'chi^2 stat: {chi2_statistic}, p_value: {p_value}')


def calc_fisher_freq_rare_by_smoking(data: pd.DataFrame) -> None:
    known_smokers = data[data[FLD_SMOKING].str.len() > 0]

    not_smokers = known_smokers[FLD_SMOKING].str.contains('не ')
    smokers = ~not_smokers

    has_mutations_condition = ~known_smokers[FLD_EGFR_TYPE].str.contains('WT')

    is_ex19del = known_smokers[FLD_EGFR_TYPE].str.contains('ex19del')
    is_l858r = known_smokers[FLD_EGFR_TYPE].str.contains('L858R')

    is_s768i = known_smokers[FLD_EGFR_TYPE].str.contains('S768I')
    is_l703v = known_smokers[FLD_EGFR_TYPE].str.contains('L703V')
    is_g779c = known_smokers[FLD_EGFR_TYPE].str.contains('G779C')
    is_d761y = known_smokers[FLD_EGFR_TYPE].str.contains('D761Y')

    freq_mutations_condition = (is_ex19del | is_l858r) & \
                               ~(is_s768i | is_l703v | is_g779c | is_d761y)

    is_double_mut = data[FLD_EGFR_TYPE].str.contains(re.escape('+'))

    has_freq_mutations = freq_mutations_condition
    has_rare_mutations = has_mutations_condition & ~freq_mutations_condition
    has_rare_double_mutations = has_rare_mutations & is_double_mut

    smokers_have_freq = known_smokers[smokers & has_freq_mutations].shape[0]
    smokers_have_rare = known_smokers[smokers & has_rare_mutations].shape[0]
    smokers_have_double_rare = known_smokers[smokers & has_rare_double_mutations].shape[0]
    not_smokers_freq = known_smokers[not_smokers & has_freq_mutations].shape[0]
    not_smokers_rare = known_smokers[not_smokers & has_rare_mutations].shape[0]
    not_smokers_double_rare = known_smokers[not_smokers & has_rare_double_mutations].shape[0]

    print('\t\tРедкие мутации\tЧастые мутации')
    print(f'Курят\t\t{smokers_have_rare}\t\t\t{smokers_have_freq}')
    print(f'Не курят\t{not_smokers_rare}\t\t\t{not_smokers_freq}')

    table = [
        [smokers_have_rare, smokers_have_freq],
        [not_smokers_rare, not_smokers_freq],
    ]

    odds_ratio, p_value = fisher_exact(table)
    print(f'Odds ratio: {odds_ratio}, p_value: {p_value}')

    chi2_statistic, p_value, dof, expected = chi2_contingency(table)
    print(f'chi^2 stat: {chi2_statistic}, p_value: {p_value}')

    # Freq & double rare

    print('\t\tЧастые мутации\tДвойные редкие мутации')
    print(f'Курят\t\t{smokers_have_freq}\t\t\t{smokers_have_double_rare}')
    print(f'Не курят\t{not_smokers_freq}\t\t\t{not_smokers_double_rare}')

    table = [
        [smokers_have_freq, smokers_have_double_rare],
        [not_smokers_freq, not_smokers_double_rare],
    ]

    odds_ratio, p_value = fisher_exact(table)
    print(f'Odds ratio: {odds_ratio}, p_value: {p_value}')

    chi2_statistic, p_value, dof, expected = chi2_contingency(table)
    print(f'chi^2 stat: {chi2_statistic}, p_value: {p_value}')


def main() -> None:
    csv = read_csv(CSV_PATH)
    compare_with_and_without_mutations(csv)
    compare_freq_and_rare_mutations(csv)
    calc_fisher_by_sex(csv)
    calc_fisher_by_smoking(csv)
    calc_fisher_freq_rare_by_sex(csv)
    calc_fisher_freq_rare_by_smoking(csv)


if __name__ == '__main__':
    main()
