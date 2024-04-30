import os
import re
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, chi2_contingency, fisher_exact, shapiro
from collections import namedtuple

SEPARATOR = ';'
CSV_PATH = 'data/oncolab_freq/list_2.csv'


FLD_EGFR_TYPE = 'EGFR type'
FLD_AGE = 'Возр'
FLD_SEX = 'Пол'
FLD_SMOKING = 'Статус курения'


Conditions = namedtuple('Conditions', ['name', 'condition'])


def read_csv(filepath: str) -> pd.DataFrame:
    data = pd.read_csv(filepath, sep=SEPARATOR)
    return data[data[FLD_AGE] > 1]


def log_to_file(data: pd.DataFrame, sheet_name: str) -> None:
    OUT_DIR = 'логи'

    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)
    with pd.ExcelWriter(f'{OUT_DIR}/{sheet_name}.xlsx', engine='openpyxl', mode='w') as writer:
        data.to_excel(writer, index=False)


def is_ex19del(data: pd.DataFrame):
    return data[FLD_EGFR_TYPE].str.contains('ex19del')


def is_l858r(data: pd.DataFrame):
     return data[FLD_EGFR_TYPE].str.contains('L858R')


def is_s768i(data: pd.DataFrame):
    return data[FLD_EGFR_TYPE].str.contains('S768I')


def is_l703v(data: pd.DataFrame):
    return data[FLD_EGFR_TYPE].str.contains('L703V')


def is_g779c(data: pd.DataFrame):
    return data[FLD_EGFR_TYPE].str.contains('G779C')


def is_d761y(data: pd.DataFrame):
    return data[FLD_EGFR_TYPE].str.contains('D761Y')


def is_ex20ins(data: pd.DataFrame):
    return data[FLD_EGFR_TYPE].str.contains('ex20ins')


def is_g719x(data: pd.DataFrame):
    return data[FLD_EGFR_TYPE].str.contains('G719')


def is_l861q(data: pd.DataFrame):
    return data[FLD_EGFR_TYPE].str.contains('L861Q')


def is_e709x(data: pd.DataFrame):
    return data[FLD_EGFR_TYPE].str.contains('E709')


def is_double_mutations(data: pd.DataFrame):
    return data[FLD_EGFR_TYPE].str.contains(re.escape('+'))


def has_mutations(data: pd.DataFrame):
    return ~(data[FLD_EGFR_TYPE].str.contains('WT'))


def has_no_mutations(data: pd.DataFrame):
    return ~(has_mutations(data))


def has_freq_mutations(data: pd.DataFrame):
    return (is_ex19del(data) | is_l858r(data)) & \
          ~(is_s768i(data) | is_l703v(data) | is_g779c(data) | is_d761y(data))


def has_rare_mutations(data: pd.DataFrame):
    return has_mutations(data) & ~(has_freq_mutations(data))


def has_rare_double_mutations(data: pd.DataFrame):
    return has_rare_mutations(data) & is_double_mutations(data)


def has_men(data: pd.DataFrame):
    return data[FLD_SEX].str.contains('м')


def has_women(data: pd.DataFrame):
    return ~(has_men(data))


def get_known_smokers(data: pd.DataFrame):
    return data[data[FLD_SMOKING].str.len() > 0]


def has_unknown_smokers(data: pd.DataFrame):
    return ~(data[FLD_SMOKING].str.len() > 0)


def has_no_smokers(data: pd.DataFrame):
    return data[FLD_SMOKING].str.contains('не кур', na=False)


def has_smokers(data: pd.DataFrame):
    return data[FLD_SMOKING].str.contains('кур') & ~has_no_smokers(data)


def calc_mann_whitneyu(patients, catigories):
    first_range = patients[catigories[0].condition][FLD_AGE]
    second_range = patients[catigories[1].condition][FLD_AGE]

    statistics, p_value = mannwhitneyu(first_range, second_range)
    print(f'{catigories[0].name} ({first_range.size}) и '
          f'{catigories[1].name} ({second_range.size}) '
          f'Mann-Whitney U statistic: {statistics}, P-value: {p_value}')


def calc_shapiro(patients, category):
    patients_by_age = patients[category.condition][FLD_AGE]
    statistics, p_value = shapiro(patients_by_age)

    print(f'Шапиро-Уилко: {category.name} '
          f'Test statistic: {statistics}, P-value: {p_value}')


def plot_distributions(patients, categories):
    plots = {category.name: patients[category.condition][FLD_AGE] for category in categories}

    distributions_by_age = pd.DataFrame(plots)
    ax = distributions_by_age.plot(kind='box')
    ax.grid(axis='y')


def calc_fisher_and_chi2_exact(patients, catigories, groups):
    group_1 = [
        patients[catigories[0].condition & groups[0].condition].shape[0],
        patients[catigories[1].condition & groups[0].condition].shape[0]
    ]

    group_2 = [
        patients[catigories[0].condition & groups[1].condition].shape[0],
        patients[catigories[1].condition & groups[1].condition].shape[0]
    ]

    table = [group_1, group_2]

    total_first_condition = patients[catigories[0].condition].shape[0]
    total_second_condition = patients[catigories[1].condition].shape[0]

    percentage_1st_condition_1st_group = get_percentage(group_1[0], total_first_condition)
    percentage_1st_condition_2st_group = get_percentage(group_2[0], total_first_condition)
    percentage_2st_condition_1st_group = get_percentage(group_1[1], total_second_condition)
    percentage_2st_condition_2st_group = get_percentage(group_2[1], total_second_condition)

    print(f'\t\t{catigories[0].name} | {catigories[1].name}\n'
          f'{groups[0].name}\t{group_1[0]} ({percentage_1st_condition_1st_group}) | '
          f'{group_1[1]} ({percentage_2st_condition_1st_group})\n'
          f'{groups[1].name}\t{group_2[0]} ({percentage_1st_condition_2st_group}) | '
          f'{group_2[1]} ({percentage_2st_condition_2st_group})')

    odds_ratio, p_value = fisher_exact(table)
    print(f'Фишер: Odds ratio: {odds_ratio}, p_value: {p_value}')

    chi2_statistic, p_value, dof, expected = chi2_contingency(table)
    print(f'chi^2 stat: {chi2_statistic}, p_value: {p_value}')


def get_percentage(count, total_count):
    return f'{count / total_count * 100.0:.2f}%'


def calc_cases_by_age(data: pd.DataFrame, min_age: int, max_age: int) -> None:
    is_greater_or_eq_min_age = data[FLD_AGE] >= min_age
    is_less_or_eq_max_age = data[FLD_AGE] <= max_age
    age_condition = is_greater_or_eq_min_age & is_less_or_eq_max_age

    all = data[FLD_AGE].size
    count = data[age_condition].shape[0]

    print(f'{min_age}-{max_age}:\t{count} ({get_percentage(count, all)})')


def print_fisher_men_and_women(patients, catigories):
    calc_fisher_and_chi2_exact(patients, catigories, [
        Conditions('Мужчины', has_men(patients)),
        Conditions('Женщины', has_women(patients))
    ])


def print_fisher_smokers_and_non_smokers(patients, catigories):
    calc_fisher_and_chi2_exact(patients, catigories, [
        Conditions('Курящие', has_smokers(patients)),
        Conditions('Некурящие', has_no_smokers(patients))
    ])


def print_statistics_for_catigories(patients, catigories):
    calc_mann_whitneyu(patients, catigories)
    print('')
    print_fisher_men_and_women(patients, catigories)
    print_fisher_smokers_and_non_smokers(patients, catigories)
    print('')


def print_statistics(patients):
    print_statistics_for_catigories(patients, [
        Conditions('Без мутаций', has_no_mutations(patients)),
        Conditions('С мутациями', has_mutations(patients))
    ])
    print_statistics_for_catigories(patients, [
        Conditions('Редкие мутации', has_rare_mutations(patients)),
        Conditions('Частые мутаци', has_freq_mutations(patients))
    ])
    print_statistics_for_catigories(patients, [
        Conditions('Частые мутаци', has_freq_mutations(patients)),
        Conditions('Редкие двойные мутации', has_rare_double_mutations(patients))
    ])
    print_statistics_for_catigories(patients, [
        Conditions('Частые мутаци', has_freq_mutations(patients)),
        Conditions('ex20', is_ex20ins(patients))
    ])
    print_statistics_for_catigories(patients, [
        Conditions('Частые мутаци', has_freq_mutations(patients)),
        Conditions('G719X', is_g719x(patients))
    ])
    print_statistics_for_catigories(patients, [
        Conditions('Частые мутаци', has_freq_mutations(patients)),
        Conditions('L861Q', is_l861q(patients))
    ])
    print_statistics_for_catigories(patients, [
        Conditions('Частые мутаци', has_freq_mutations(patients)),
        Conditions('S768I', is_s768i(patients))
    ])
    print_statistics_for_catigories(patients, [
        Conditions('Частые мутаци', has_freq_mutations(patients)),
        Conditions('E709X', is_e709x(patients))
    ])
    print_statistics_for_catigories(patients, [
        Conditions('Частые мутаци', has_freq_mutations(patients)),
        Conditions('ex19del', is_ex19del(patients))
    ])
    print_statistics_for_catigories(patients, [
        Conditions('Частые мутаци', has_freq_mutations(patients)),
        Conditions('L858R', is_l858r(patients))
    ])


def print_cases_rare_by_ages(title, patients):
    total_count = patients.shape[0]
    patient_ages = patients[FLD_AGE]
    print(f'{title} - Медиана: {patient_ages.mean()}, '
          f'Диапазон: ({patient_ages.min()}-{patient_ages.max()})')

    calc_cases_by_age(patients, 0, 40)
    calc_cases_by_age(patients, 41, 50)
    calc_cases_by_age(patients, 51, 60)
    calc_cases_by_age(patients, 61, 70)
    calc_cases_by_age(patients, 71, 999)

    men = patients[has_men(patients)].shape[0]
    women = patients[has_women(patients)].shape[0]
    print(f'Мужчины: {men} ({get_percentage(men, total_count)}')
    print(f'Женщины: {women} ({get_percentage(women, total_count)})')

    smokers = patients[has_smokers(patients)].shape[0]
    non_smokers = patients[has_no_smokers(patients)].shape[0]
    unknown_smokers = patients[has_unknown_smokers(patients)].shape[0]
    print(f'Курящие: {smokers} ({get_percentage(smokers, total_count)}')
    print(f'Некурящие: {non_smokers} ({get_percentage(non_smokers, total_count)})')
    print(f'Неизвестно: {unknown_smokers} ({get_percentage(unknown_smokers, total_count)})')

    print(f'Всего: {total_count}')
    print('')


def print_cases_by_ages(patients):
    print('Наиболее часто встречающиеся случаи редких мутаций гена EGFR')
    print_cases_rare_by_ages('Редкие', patients[has_rare_mutations(patients)])
    print_cases_rare_by_ages('ex20ins', patients[is_ex20ins(patients)])
    print_cases_rare_by_ages('G719X', patients[is_g719x(patients)])
    print_cases_rare_by_ages('L861Q', patients[is_l861q(patients)])
    print_cases_rare_by_ages('S768I', patients[is_s768i(patients)])
    print_cases_rare_by_ages('E709X', patients[is_e709x(patients)])
    print_cases_rare_by_ages('ex19del', patients[is_ex19del(patients)])
    print_cases_rare_by_ages('L858R', patients[is_l858r(patients)])
    print('')


def print_shapiro(patients):
    calc_shapiro(patients, Conditions('Есть мутации', has_mutations(patients)))
    calc_shapiro(patients, Conditions('Без мутаций', has_no_mutations(patients)))
    calc_shapiro(patients, Conditions('Частые мутации', has_freq_mutations(patients)))
    calc_shapiro(patients, Conditions('Редкие мутации', has_rare_mutations(patients)))

    calc_shapiro(patients, Conditions('Редкие двойные мутации', has_rare_double_mutations(patients)))
    calc_shapiro(patients, Conditions('ex20', is_ex20ins(patients)))
    calc_shapiro(patients, Conditions('G719X', is_g719x(patients)))
    calc_shapiro(patients, Conditions('L861Q', is_l861q(patients)))
    calc_shapiro(patients, Conditions('S768I', is_s768i(patients)))
    calc_shapiro(patients, Conditions('E709X', is_e709x(patients)))


def show_plots(patients):
    plot_distributions(patients, [
        Conditions('Есть мутации', has_mutations(patients)),
        Conditions('Без мутаций', has_no_mutations(patients))
    ])
    plot_distributions(patients, [
        Conditions('Частые мутации', has_freq_mutations(patients)),
        Conditions('Редкие мутации', has_rare_mutations(patients)),
    ])
    plot_distributions(patients, [
        Conditions('Частые мутации', has_freq_mutations(patients)),
        Conditions('Редкие двойные мутации', has_rare_double_mutations(patients)),
    ])
    plot_distributions(patients, [
        Conditions('Частые мутации', has_freq_mutations(patients)),
        Conditions('ex20', is_ex20ins(patients)),
        Conditions('G719X', is_g719x(patients)),
        Conditions('L861Q', is_l861q(patients)),
        Conditions('S768I', is_s768i(patients)),
        Conditions('E709X', is_e709x(patients)),
    ])

    plt.show()


def main() -> None:
    patients = read_csv(CSV_PATH)

    print_statistics(patients)
    # print_cases_by_ages(patients)
    print_shapiro(patients)
    show_plots(patients)
    plt.show()


if __name__ == '__main__':
    main()
