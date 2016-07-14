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


# These will be removed after all have been converted to .csv files.
# TODO (Figure out what some of these are...)
SOURCE_PATH_MAT = os.path.join('golden', 'P_ritter.mat')
DATA_PATH_MAT = os.path.join('golden', 'data_all_mean.mat')
BODYWEIGHT_PATH_MAT = os.path.join('golden', 'bw_cz.mat')  # bodyweight data
LIPID_FRACTION_PATH_MAT = os.path.join('golden', 'flip_cz.mat')

BODYWEIGHT_AND_LIPID_FRACTION_ = os.path.join(
    'monthly_bodyweight_and_lipid_fraction.csv')
BIOMONITORING_DATA_ELIMINATION_ = os.path.join('mikes_et_al_2012.csv')
BIOMONITORING_DATA_ELIMINATION_UNITS = ['ng/kg_bw/d', 'years']  # just to keep track
_CONGENER_START_PEAK_AGE_GROUP = os.path.join('cogener_start_peak_age_group.csv')
_CONGENERS_TO_EVALUATE = ['hcb', 'ppddt', 'ppdde', 'betahch', 'gammahch', 'pcb138']
_FITTING_APPROACH = ['intake', 'intake_k_elim', 'original', 't_elim']


N_MALES_ = 100
N_FEMALES_ = 100
STUDY_START_YEAR_ = 1994
STUDY_END_YEAR_ = 2009
AGE_MAX_IN_YEARS_ = 80  # [years]
ABSORP_FACTOR_ = 0.9  # absorbption factor [-]
AVERAGE_LACT_TIME_MONTHS_ = 6  # months breastfeeding (generic source?)
AGE_MOM_IN_YEARS_ = 25
AGE_GROUPINGS_ = range(15, 45, 5)  # start, stop, step
PEAK_DOSING_INTENSITY_ = 80  # [ng/kg/d]
N_OPTIMIZATION_RUNS_ = 1
REF_FEMALE_BODYWEIGHT_KG_ = 69
KINETICS_ORDER_ = 1
AVERAGE_BIRTH_LIPID_FRACTION_ = 0.215
AVERAGE_BIRTH_BODYWEIGHT_KG_ = 3.1744


def main(n_optimization_runs=N_OPTIMIZATION_RUNS_,
         biomonitoring_data=BIOMONITORING_DATA_ELIMINATION_,
         congener_start_peak_age_group=_CONGENER_START_PEAK_AGE_GROUP,
         congeners_to_evaluate=_CONGENERS_TO_EVALUATE,  # list
         kinetics_order=KINETICS_ORDER_,
         plot_kinetics=False,
         n_male=N_MALES_,
         n_female=N_FEMALES_,
         study_start_year=STUDY_START_YEAR_,
         study_end_year=STUDY_END_YEAR_,
         average_lact_time_months_=AVERAGE_LACT_TIME_MONTHS_,
         age_max_in_years=AGE_MAX_IN_YEARS_,
         age_mom_in_years=AGE_MOM_IN_YEARS_,
         absorbption_factor=ABSORP_FACTOR_,
         average_birth_lipid_fraction_=AVERAGE_BIRTH_LIPID_FRACTION_,
         average_birth_bodyweight_kg_=AVERAGE_BIRTH_BODYWEIGHT_KG_,
         bodyweight_and_lipid_fraction=BODYWEIGHT_AND_LIPID_FRACTION_):

    # error handling calls
    def biomonitoring_data_error_handling(a, b):
        eval_bio_diff = list(set(a) - set(b))

        if not eval_bio_diff:
            print '''The following chemicals with biomonitoring data will not be evaluated'''.replace(
                '/t', ''), eval_bio_diff
        else:
            print ''' Biomonitoring data may be available, but not evaluating chemicals:'''.replace(
                '/t', ''), eval_bio_diff

    # kinetics - import biomonitoring data if needed
    # TODO(build in argparse here for csv extension)
    kinetics_from_biomonitoring_data = []
    if biomonitoring_data:
        biomonitoring_data = pd.read_csv(biomonitoring_data)
        biomonitoring_data_time = biomonitoring_data.ix[:, 0]
        congeners_to_evaluate = congeners_to_evaluate
        biomonitoring_data_congoners = list(biomonitoring_data[-1:])

        biomonitoring_data_error_handling(
            congeners_to_evaluate, biomonitoring_data_congoners)

        colors = cm.rainbow(np.linspace(0, 1, len(congeners_to_evaluate)))
        for congener, c in zip(congeners_to_evaluate, colors):
            if kinetics_order == 1:
                x = np.array(biomonitoring_data_time)
                y = np.log10(np.array(biomonitoring_data[congener]))
                (slope, intercept, r_value, p_value, std_err) = stats.linregress(x, y)

                # TODO(output the kinetic rates and fit values to a table)
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
                y_regression_line = polyval([slope, intercept], x)
                kinetics_from_biomonitoring_data.append(slope)
                if plot_kinetics:

                    plt.plot(10**x, 10**y, color=c,
                             marker='o', linestyle='', label=str(
                                 'biomonitoring data; ' + congener))
                    plt.plot(10**x, 10**y_regression_line,
                             color=c, marker='', linestyle='-')

        if plot_kinetics:
            plt.show()

    # assigning chemical specific parameters
    congener_start_peak_age_group_df = pd.read_csv(congener_start_peak_age_group)
    biomonitoring_data_error_handling(congeners_to_evaluate,
                                      congener_start_peak_age_group_df)

    # only use values which are in the target cogeners
    congener_start_peak_age_group_df[congener_start_peak_age_group_df.isin(
        congeners_to_evaluate)]

    # calculate mass of lipids for each month
    def calculate_mass_lipids(row):
        mass_lipids = row.avg_monthly_bodyweight_kg * row.lipid_fraction
        return mass_lipids

    def calculate_kinetic_growth_from_lipid_frac(row):
        kinetic_growth = (
            row.mass_lipids_kg - row.lipid_frac_offset) / row.mass_lipids_kg
        return kinetic_growth

    if bodyweight_and_lipid_fraction:
        bodyweight_and_lipid_fraction = pd.read_csv(bodyweight_and_lipid_fraction)
        bodyweight_and_lipid_fraction[
            'mass_lipids_kg'] = bodyweight_and_lipid_fraction.apply(calculate_mass_lipids, 1)

        bodyweight_and_lipid_fraction['lipid_frac_offset'] = bodyweight_and_lipid_fraction[
            'mass_lipids_kg'].shift(1)
        bodyweight_and_lipid_fraction['lipid_frac_offset'].iloc[
            0] = average_birth_lipid_fraction_ * average_birth_bodyweight_kg_
        bodyweight_and_lipid_fraction['k_growth'] = bodyweight_and_lipid_fraction.apply(
            calculate_kinetic_growth_from_lipid_frac, 1)
    else:
        print("please include .csv file containing the'average_body_weight_in_kg' in the first column and the 'lipid_fraction' in the second column")

    # calculating breastfeeding kinetics

    # [months] age of mother at birth in months
    age_mom_in_months = age_mom_in_years * 12

    # [years] age in years with interval of 1 year
    age_in_years_array = range(1, age_max_in_years)  # start, stop, step=1

    # [months] age in weeks with interval of 1 month
    age_m = range(1, max_age_in_months)  # start, stop, step=1

    # [years] age in years with interval of 1 month
    age_years_in_months = np.linspace(start=1.0 / 12,
                                      stop=age_max_in_years,
                                      num=age_max_in_years * 12,
                                      endpoint=True)
    lact_kinetics = pd.DataFrame()
    lact_kinetics['age_index_years'] = age_years_in_months
    lact_kinetics['age_index_months'] = lact_kinetics['age_index_years'] * 12

    print lact_kinetics.head()
    p_source = mat_to_df(SOURCE_PATH_MAT).T
    all_mean_data = mat_to_df(DATA_PATH_MAT).to_csv('test.csv')


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
    # age iteration parameters
    # Time variables (single-value parameters --> input_single matrix)
    years_to_months = 12  # [months/year] unit conversion factor
    month_to_days = 30  # for loading characterization

    # TODO (eli) calendarize this function?
    months_to_days = 30
    n_people = n_male + n_female


    # TODO(peak dosing intensity on/off? set function?)
    peak_dosing_intensity = 80 * absorbption_factor / 1000 * month_to_days

    # print age_m

    # set the kinetics for females that give birth
    # for female in

    # mat_to_df(BODYWEIGHT_PATH_MAT).to_csv(
    #     'avg_monthly_bodyweight_kg.csv')
    # bodyweight_data = pd.read_csv('avg_monthly_bodyweight_kg.csv')
    # bodyweight_data_2 = pd.DataFrame()
    # bodyweight_data_2['avg_monthly_bodyweight_kg'] = bodyweight_data['0']
    # print bodyweight_data_2
    # bodyweight_data_2['avg_monthly_bodyweight_kg'].to_csv(
    #     'avg_monthly_bodyweight_kg.csv', index=False, header=True)

    # mat_to_df(LIPID_FRACTION_PATH_MAT).to_csv(
    #     'lipid_fraction.csv')
    # lipid_fraction = pd.read_csv('lipid_fraction.csv')
    # lipid_fraction_2 = pd.DataFrame()
    # lipid_fraction_2['lipid_fraction'] = lipid_fraction['0'] / 100.
    # print lipid_fraction_2
    # lipid_fraction_2['lipid_fraction'].to_csv(
    #     'lipid_fraction.csv', index=False, header=True)

    # b = pd.read_csv('lipid_fraction.csv')
    # a = pd.read_csv('avg_monthly_bodyweight_kg.csv')
    # monthly_bodyweight_and_lipid_fraction = pd.concat([a, b], 1)
    # monthly_bodyweight_and_lipid_fraction.to_csv(
    #     'monthly_bodyweight_and_lipid_fraction.csv', index=False, header=True)

"""
