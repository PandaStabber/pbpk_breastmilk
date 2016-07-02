# coding=utf-8
import matplotlib.pyplot as plt
import scipy.io
import pandas as pd
from scipy.optimize import curve_fit
from scipy import polyval, stats
import numpy as np
import os
import matplotlib.cm as cm


def mat_to_df(path):
    mat = scipy.io.loadmat(path,
                           matlab_compatible=True,
                           variable_names=None)
    mat_heads_init = list(mat)[0]
    mat_matrix = mat[mat_heads_init]
    df = pd.DataFrame(mat_matrix)
    return df


# default source files
SOURCE_PATH_MAT = os.path.join('golden', 'P_ritter.mat')

# default data
DATA_PATH_MAT = os.path.join('golden', 'data_all_mean.mat')

# Chem elimination data (ng/glip)
# organized: time, conjoiner1, conjoiner2 (ng/kgbw/d), ...
BIOMONITORING_DATA_ELIMINATION = os.path.join('mikes_et_al_2012.csv')
BIOMONITORING_DATA_ELIMINATION_UNITS = ['ng/kg_bw/d', 'years']

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
REF_FEMALE_BODY_WEIGHT_KG_ = 69
KINETICS_ORDER_ = 1


_CONGENER_START_PEAK_AGE_GROUP = os.path.join('cogener_start_peak_age_group.csv')

_CONGENERS_TO_EVALUATE = ['hcb', 'ppddt', 'ppdde', 'betahch', 'gammahch', 'pcb138']
_FITTING_APPROACH = ['intake', 'intake_k_elim', 'original', 't_elim']


def congener_SPA_lookup(row):
    return _CONGENER_START_PEAK_AGE_GROUP.get(row.congener_SPA_lookup,
                                              [0, 0])


def main(n_optimization_runs=N_OPTIMIZATION_RUNS_,
         biomonitoring_data=BIOMONITORING_DATA_ELIMINATION,
         biomonitoring_data_elimination_units=BIOMONITORING_DATA_ELIMINATION_UNITS,
         congener_start_peak_age_group=_CONGENER_START_PEAK_AGE_GROUP,
         congeners_to_evaluate=_CONGENERS_TO_EVALUATE,  # list
         kinetics_order=KINETICS_ORDER_,
         plot_kinetics=False,
         n_male=N_MALES_,
         n_female=N_FEMALES_,
         study_start_year=STUDY_START_YEAR_,
         study_end_year=STUDY_END_YEAR_,
         t_lact=T_LACT_,
         age_max_in_years=AGE_MAX_IN_YEARS_,
         age_mom_in_years=AGE_MOM_IN_YEARS_,
         absorbption_factor=ABSORP_FACTOR_):

    # kinetics - import biomonitoring data if needed
    # TODO(build in argparse here for csv extension)
    # kinetics_from_biomonitoring_data = pd.DataFrame()
    kinetics_from_biomonitoring_data = []
    if biomonitoring_data:
        biomonitoring_data = pd.read_csv(biomonitoring_data)
        biomonitoring_data_time = biomonitoring_data.ix[:, 0]
        congeners_to_evaluate = congeners_to_evaluate
        biomonitoring_data_congoners = list(biomonitoring_data[-1:])
        eval_bio_diff = list(
            set(congeners_to_evaluate) -
            set(biomonitoring_data_congoners))
        if not eval_bio_diff:
            print "All chemicals with biomonitoring data will be evaluated"
        else:
            print "biomonitoring data available, but not evaluating chemicals: ", eval_bio_diff

        colors = cm.rainbow(np.linspace(0, 1, len(congeners_to_evaluate)))
        for congener, c in zip(congeners_to_evaluate, colors):
            if kinetics_order == 1:
                x = np.array(biomonitoring_data_time)
                y = np.log10(np.array(biomonitoring_data[congener]))
                (slope, intercept, r_value, p_value, std_err) = stats.linregress(x, y)
                # TODO(output the kinetic rates and fit values to a table)
                # print(congener, " has a 1st order kinetic elimination estimate of ",
                #       slope, str('per ' + BIOMONITORING_DATA_ELIMINATION_UNITS[1]))
                # print(congener, " has an r^2 fit value of ", r_value**2)

                y_regression_line = polyval([slope, intercept], x)
                kinetics_from_biomonitoring_data.append(slope)
                if plot_kinetics:

                    plt.plot(x, 10**y, color=c,
                             marker='o', linestyle='', label=str(
                                 'biomonitoring data; ' + congener))
                    plt.plot(x, 10**y_regression_line,
                             color=c, marker='', linestyle='-')

            elif kinetics_order == 2:
                x = np.log10(np.array(biomonitoring_data_time))
                y = np.log10(np.array(biomonitoring_data[congener]))
                (slope, intercept, r_value, p_value, std_err) = stats.linregress(x, y)
                # TODO(output the kinetic rates and fit values to a table)
                # print(congener, " has a 1st order kinetic elimination estimate of ",
                #       slope, str('per ' + BIOMONITORING_DATA_ELIMINATION_UNITS[1]))
                # print(congener, " has an r^2 fit value of ", r_value**2)

                y_regression_line = polyval([slope, intercept], x)
                kinetics_from_biomonitoring_data.append(slope)
                if plot_kinetics:

                    plt.plot(10**x, 10**y, color=c,
                             marker='o', linestyle='', label=str(
                                 'biomonitoring data; ' + congener))
                    plt.plot(10**x, 10**y_regression_line,
                             color=c, marker='', linestyle='-')

        print kinetics_from_biomonitoring_data
        if plot_kinetics:
            plt.show()

    # assigning chemical specific parameters
    # import start peak age group data
    congener_start_peak_age_group_df = pd.read_csv(congener_start_peak_age_group)
    print congener_start_peak_age_group_df
    congener_start_peak_age_group_df = congener_start_peak_age_group_df[
        congener_start_peak_age_group_df.cogener in congeners_to_evaluate]

    p_source = mat_to_df(SOURCE_PATH_MAT).T
    all_mean_data = mat_to_df(DATA_PATH_MAT).to_csv('test.csv')

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

    # TODO(peak dosing intensity on/off? set function?)
    peak_dosing_intensity = 80 * absorbption_factor / 1000 * month_to_days

    # print age_m

    # set the kinetics for females that give birth
    # for female in


main()


"""[summary]
# TODO(Create function to sum up exposure from food)
As is, it's not clear how the exposure is calculated from data
D = C*IR*AF*EF*CF/BW
D = Exposure dose
C = contaminant concentration (mg/kg)
IR = intake rate of contaminated medium (mg/d)
AF = bioavailability factor (assume 1)
EF = exposure factor (assume 1)
CF = conversion factor (10e-6 kg/mg)
BW = body weight (kg)

-- EF = (F*ED)/AT, where F = frequency of exposure (d/year)
                  ED = exposure duration (years)
                  AT = averaging time = (ED*365d/yeaer)
Ref: https://www.atsdr.cdc.gov/hac/phamanual/appg.html

[description]
"""

"""[summary]
    for congener in congeners_to_evaluate:
        ss_val = congener_start_peak_age_group.get(congener)
        print ss_val

    for congeners in congener_start_peak_age_group:
        if congeners in congeners_to_evaluate:
            ss_val = congener_start_peak_age_group.get(congeners)
            print ss_val
            for item in ss_val:
                congener_start_peak_age_group_df[str(item)] = item

    print congener_start_peak_age_group_df
[description]
"""
