import scipy.io
import pandas as pd
import numpy as np


def mat_to_df(path):

    mat = scipy.io.loadmat(path,
                           matlab_compatible=True,
                           variable_names=None)
    mat_heads_init = list(mat)[0]
    mat_matrix = mat[mat_heads_init]
    df = pd.DataFrame(mat_matrix)
    return df


p_source = mat_to_df("golden_files/P_ritter.mat").T

all_mean_data = mat_to_df("golden_files/data_all_mean.mat")

N_MALES_ = 100
N_FEMALES_ = 100
STUDY_START_YEAR_ = 1994
STUDY_END_YEAR_ = 2009
AGE_MAX_IN_YEARS_ = 80  # [years]
ABSORP_FACTOR_ = 0.9  # absorbption factor [-]
T_LACT_ = 0.5  # years breastfeeding (generic source?)
AGE_MOM_IN_YEARS_ = 25
AGE_GROUPINGS_ = range(15, 45, 5)  # start, stop, step
PEAK_DOSING_INTENSITY_ = 80  # [ng/kg/d]
N_OPTIMIZATION_RUNS_ = 1

# congener organization
# start year, peak year
# end year is parameterized as a functinon of # of people
# TODO(eli) Tenzing has a time array as well, remove
_CONGENER_START_PEAK_AGE_GROUP = {
    'bHCH': [1921, 1977, 7],
    'DDE': [1921, 1974, 9],
    'DDT': [1921, 1977, 10],
    'gHCH': [1921, 1977, 7],
    'HCB': [1921, 1977, 7],
    'PCB28': [1921, 1975, 7],
    'PCB52': [1921, 1977, 7],
    'PCB101': [1921, 1977, 7],
    'PCB118': [1921, 1977, 7],
    'PCB138': [1921, 1977, 7],
    'PCB153': [1921, 1977, 7],
    'PCB170': [1921, 1977, 7],
    'PCB180': [1921, 1977, 6]}

_CONGENERS_TO_EVALUATE = ['PCB180', 'bHCH', 'DDE', 'DDT']


def congener_start_peak_lookup(row):
    """Find valence.
    notes: match ENM and return valence
    :param row:
    :return: corresponding electrolyte valence
    """
    return _CONGENER_START_AND_PEAK.get(row.congener_start_peak_lookup, [0, 0])


def congener_age_group_lookup(row):
    return _CONGENER_AGE_GROUP_DEFINITION.get(
        row.congener_age_group_lookup, [0, 0, 0])


def main(n_optimization_runs=N_OPTIMIZATION_RUNS_,
         n_male=N_MALES_,
         n_female=N_FEMALES_,
         study_start_year=STUDY_START_YEAR_,
         study_end_year=STUDY_END_YEAR_,
         t_lact=T_LACT_,
         age_max_in_years=AGE_MAX_IN_YEARS_,
         age_mom_in_years=AGE_MOM_IN_YEARS_,
         absorbption_factor=ABSORP_FACTOR_,
         congener_start_peak_age_group=_CONGENER_START_PEAK_AGE_GROUP,
         congeners_to_evaluate=_CONGENERS_TO_EVALUATE):

    # age iteration parameters
    # Time variables (single-value parameters --> input_single matrix)
    years_to_months = 12  # [months/year] unit conversion factor
    month_to_days = 30  # for loading characterization

    # TODO (eli) calendarize this function?
    months_to_days = 30
    n_people = n_male + n_female

    # [months] maximum age of a person in months
    max_age_in_months = age_max_in_years * years_to_months

    # [months] age of mother at birth in months
    age_mom_in_months = age_mom_in_years * years_to_months

    # [years] age in years with interval of 1 year
    age_in_years_array = range(1, age_max_in_years)  # start, stop, step=1

    # [months] age in weeks with interval of 1 month
    age_m = range(1, max_age_in_months)  # start, stop, step=1

    # [years] age in years with interval of 1 month
    age_years_in_months = np.linspace(start=1.0 / years_to_months,
                                      stop=age_max_in_years,
                                      num=age_max_in_years * years_to_months,
                                      endpoint=True)

    # TODO(eli) peak dosing intensity on/off? set function?
    # TODO(eli) ensure 80 is ng/kg/d
    peak_dosing_intensity = 80 * absorbption_factor / 1000 * month_to_days

    congener_start_peak_age_group_df = pd.DataFrame(
        congeners_to_evaluate, columns=['congener_name'])

    for congeners in congener_start_peak_age_group:
        if congeners in congeners_to_evaluate:
            ss_val = congener_start_peak_age_group.get(congeners)
            for item in ss_val:
                congener_start_peak_age_group_df[str(item)] = item

    print congener_start_peak_age_group_df

main()
