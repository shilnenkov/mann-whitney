import pandas as pd
import os
import re
from scipy.stats import mannwhitneyu, chi2_contingency, fisher_exact
from collections import namedtuple

SEPARATOR = ';'
CSV_PATH = 'data/oncolab_freq/list_2.csv'


FLD_EGFR_TYPE = 'EGFR type'
FLD_AGE = 'Возр'
FLD_SEX = 'Пол'
FLD_SMOKING = 'Статус курения'


FisherConditions = namedtuple('FisherConditions', ['name', 'condition'])

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


def has_no_smokers(data: pd.DataFrame):
    return data[FLD_SMOKING].str.contains('не ')


def has_smokers(data: pd.DataFrame):
    return ~(has_no_smokers(data))


def calc_mann_whitneyu(first_range_name, first_range, second_range_name, second_range):
    statistics, p_value = mannwhitneyu(first_range, second_range)
    print(f'{first_range_name} ({first_range.size}) и '
          f'{second_range_name} ({second_range.size}) '
          f'Mann-Whitney U statistic: {statistics}, P-value: {p_value}')


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

    print(f'\t\t{catigories[0].name}\t{catigories[1].name}\n'
          f'{groups[0].name}\t{group_1[0]} | {group_1[1]}\n'
          f'{groups[1].name}\t{group_2[0]} | {group_2[1]}')

    odds_ratio, p_value = fisher_exact(table)
    print(f'Фишер: Odds ratio: {odds_ratio}, p_value: {p_value}')

    chi2_statistic, p_value, dof, expected = chi2_contingency(table)
    print(f'chi^2 stat: {chi2_statistic}, p_value: {p_value}')


def get_percentage(count, total_count):
    return f'{count / total_count * 100.0:.2f}%'


def calc_cases_by_age(data: pd.DataFrame, min_age: int, max_age: int) -> None:
    is_greater_or_eq_min_age = data[FLD_AGE] > min_age
    is_less_or_eq_max_age = data[FLD_AGE] <= max_age
    age_condition = is_greater_or_eq_min_age & is_less_or_eq_max_age

    all = data[FLD_AGE].size
    count = data[age_condition].shape[0]

    print(f'{min_age}-{max_age}:\t{count} ({get_percentage(count, all)})')


def print_with_and_without_mutations(patients):
    known_smokers = get_known_smokers(patients)

    calc_mann_whitneyu(
        'С мутациями', patients[has_mutations(patients)][FLD_AGE],
        'Без мутаций', patients[has_no_mutations(patients)][FLD_AGE]
    )
    print('')
    calc_fisher_and_chi2_exact(patients,
                               [
                                   FisherConditions('Без мутаций', has_no_mutations(patients)),
                                   FisherConditions('С мутациями', has_mutations(patients))
                               ],
                               [
                                   FisherConditions('Мужчины', has_men(patients)),
                                   FisherConditions('Женщины', has_women(patients))
                               ]
                               )
    calc_fisher_and_chi2_exact(known_smokers,
                               [
                                   FisherConditions('Без мутаций', has_no_mutations(known_smokers)),
                                   FisherConditions('С мутациями', has_mutations(known_smokers))
                               ],
                               [
                                   FisherConditions('Курят', has_smokers(known_smokers)),
                                   FisherConditions('Не курят', has_no_smokers(known_smokers))
                               ]
                               )
    print('')


def print_freq_and_rare_mutations(patients):
    known_smokers = get_known_smokers(patients)

    calc_mann_whitneyu(
        'Частые', patients[has_freq_mutations(patients)][FLD_AGE],
        'редкие', patients[has_rare_mutations(patients)][FLD_AGE]
    )
    print('')
    calc_fisher_and_chi2_exact(patients,
                               [
                                   FisherConditions('Редкие мутации', has_rare_mutations(patients)),
                                   FisherConditions('Частые мутаци', has_freq_mutations(patients))
                               ],
                               [
                                   FisherConditions('Мужчины', has_men(patients)),
                                   FisherConditions('Женщины', has_women(patients))
                               ]
                               )

    calc_fisher_and_chi2_exact(known_smokers,
                               [
                                   FisherConditions('Редкие мутации', has_rare_mutations(known_smokers)),
                                   FisherConditions('Частые мутаци', has_freq_mutations(known_smokers))
                               ],
                               [
                                   FisherConditions('Курят', has_smokers(known_smokers)),
                                   FisherConditions('Не курят', has_no_smokers(known_smokers))
                               ]
                               )
    print('')


def print_freq_and_rare_double_mutations(patients):
    known_smokers = get_known_smokers(patients)

    calc_mann_whitneyu(
        'Частые', patients[has_freq_mutations(patients)][FLD_AGE],
        'редкие двойные', patients[has_rare_double_mutations(patients)][FLD_AGE]
    )
    print('')
    calc_fisher_and_chi2_exact(patients,
                               [
                                   FisherConditions('Частые мутации', has_freq_mutations(patients)),
                                   FisherConditions('Редкие двойные мутации', has_rare_double_mutations(patients))
                               ],
                               [
                                   FisherConditions('Мужчины', has_men(patients)),
                                   FisherConditions('Женщины', has_women(patients))
                               ]
                               )

    calc_fisher_and_chi2_exact(known_smokers,
                               [
                                   FisherConditions('Частые мутации', has_freq_mutations(known_smokers)),
                                   FisherConditions('Редкие двойные мутации', has_rare_double_mutations(known_smokers))
                               ],
                               [
                                   FisherConditions('Курят', has_smokers(known_smokers)),
                                   FisherConditions('Не курят', has_no_smokers(known_smokers))
                               ]
                               )
    print('')


def print_cases_rare_by_ages(title, patients):
    total_count = patients.shape[0]
    print(f'{title} - Медиана: {patients[FLD_AGE].mean()}')

    calc_cases_by_age(patients, 0, 40)
    calc_cases_by_age(patients, 41, 50)
    calc_cases_by_age(patients, 51, 60)
    calc_cases_by_age(patients, 61, 70)
    calc_cases_by_age(patients, 70, 999)

    men = patients[has_men(patients)].shape[0]
    women = patients[has_women(patients)].shape[0]
    print(f'Мужчины: {men} ({get_percentage(men, total_count)}')
    print(f'Женщины: {women} ({get_percentage(women, total_count)})')

    known_smokers = get_known_smokers(patients)
    smokers = known_smokers[has_smokers(known_smokers)].shape[0]
    non_smokers = known_smokers[has_no_smokers(known_smokers)].shape[0]
    unknown_smokers = total_count - known_smokers.shape[0]
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


def main() -> None:
    patients = read_csv(CSV_PATH)

    print_with_and_without_mutations(patients)
    print_freq_and_rare_mutations(patients)
    print_freq_and_rare_double_mutations(patients)
    print_cases_by_ages(patients)


if __name__ == '__main__':
    main()
