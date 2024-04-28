import pandas as pd
from scipy.stats import mannwhitneyu

SEPARATOR = ';'
CSV_PATH = 'data/oncolab_freq/list_2.csv'

FLD_EGFR_TYPE = 'EGFR type'
FLD_AGE = 'Возр'

MUTATATION_WT = 'WT'
MUTATATION_EX19DEL = 'ex19del'
MUTATATION_L858R = 'L858R'


def read_csv(filepath: str) -> pd.DataFrame:
    return pd.read_csv(filepath, sep=SEPARATOR)


def log_to_file(data: pd.DataFrame, filename: str) -> None:
    data.to_csv(filename, index=False)


def compare_with_and_without_mutations(data: pd.DataFrame) -> None:
    has_no_mutations_condition = data[FLD_EGFR_TYPE].str.contains(MUTATATION_WT)
    has_mutations_condition = ~(has_no_mutations_condition)

    with_mutations = data[has_mutations_condition][FLD_AGE]
    without_mutations = data[has_no_mutations_condition][FLD_AGE]

    statistics, p_value = mannwhitneyu(with_mutations, without_mutations)
    print(f'С мутациями ({with_mutations.size}) и '
          f'без ({without_mutations.size}) - '
          f'Mann-Whitney U statistic: {statistics}, P-value: {p_value}')

    log_to_file(data[has_mutations_condition], "with_mutations.csv")
    log_to_file(data[has_no_mutations_condition], "without_mutations.csv")


def compare_freq_and_rare_mutations(data: pd.DataFrame) -> None:
    has_mutations_condition = ~(data[FLD_EGFR_TYPE].str.contains(MUTATATION_WT))

    is_ex19del = data[FLD_EGFR_TYPE].str.contains(MUTATATION_EX19DEL)
    is_l858r = data[FLD_EGFR_TYPE].str.contains(MUTATATION_L858R)
    freq_mutations_condition = is_ex19del | is_l858r

    has_freq_mutations = freq_mutations_condition
    has_rare_mutations = has_mutations_condition & ~(freq_mutations_condition)

    freq_mutations = data[has_freq_mutations][FLD_AGE]
    rare_mutations = data[has_rare_mutations][FLD_AGE]

    statistics, p_value = mannwhitneyu(freq_mutations, rare_mutations)
    print(f'Частые ({freq_mutations.size}) и '
          f'редкие ({rare_mutations.size}) - '
          f'Mann-Whitney U statistic: {statistics}, P-value: {p_value}')

    log_to_file(data[has_freq_mutations], "freq_mutations.csv")
    log_to_file(data[has_rare_mutations], "rare_mutations.csv")


def main() -> None:
    csv = read_csv(CSV_PATH)
    compare_with_and_without_mutations(csv)
    compare_freq_and_rare_mutations(csv)


if __name__ == '__main__':
    main()
