# coding=utf-8
import os
from functions import *

SOURCE_PATH_MAT = os.path.join('golden', 'P_ritter.mat')
DATA_PATH_MAT = os.path.join('golden', 'data_all_mean.mat')
BODYWEIGHT_PATH_MAT = os.path.join('golden', 'bw_cz.mat')  # bodyweight data
LIPID_FRACTION_PATH_MAT = os.path.join('golden', 'flip_cz.mat')

BODYWEIGHT_AND_LIPID_FRACTION_ = os.path.join(
    'monthly_bodyweight_and_lipid_fraction.csv')
BIOMONITORING_DATA_ELIMINATION_ = os.path.join('mikes_et_al_2012.csv')
BIOMONITORING_DATA_ELIMINATION_UNITS = [
    'ng/kg_bw/d', 'years']  # just to keep track
_cong_spag = os.path.join(
    'cogener_start_peak_age_group.csv')
cte_ = ['hcb', 'ppddt',
        'ppdde', 'betahch', 'gammahch', 'pcb138']

_CONGENERS_HALFLIFE_DICT = {'hcb': [10, np.log(2) / 5, np.log(2) / 5],
                            'ppddt': [5, np.log(2) / 5, np.log(2) / 5]}

cte_ = list(_CONGENERS_HALFLIFE_DICT.keys())

TIMESTEPS_PER_MONTH = 1
POP_START_YEAR_ = 0
STUDY_START_YEAR_ = 0
STUDY_END_YEAR_ = 100
AGE_MAX_IN_YEARS_ = 80  # [years]
ABSORP_FACTOR_ = 0.9  # absorbption factor [-]
AVERAGE_LACT_TIME_ = [0.5] * 5  # months breastfeeding (generic source?)
K_LAC_GEN = [np.log(2)/5] * 5
AGE_MOM_IN_YEARS_ = 25

PEAK_DOSING_INTENSITY_ = 80  # [ng/kg/d]
REF_FEMALE_BODYWEIGHT_KG_ = 69
AVERAGE_MOTHER_LIPID_FRACTION_ = 0.35  # I MADE THIS UP - DEFAULT CHANGE
KINETICS_ORDER_ = 1
AVERAGE_BIRTH_LIPID_FRACTION_ = 0.215
body_w_ave = 3.1744
K_ELIM = np.log(2) / 5
PEAK_YEAR = [35]

def pk_model(bm_data=BIOMONITORING_DATA_ELIMINATION_,
             cong_spag=_cong_spag,
             cte=cte_,
             kinetics_order=KINETICS_ORDER_,
             gens=1,
             plot_kinetics=False,
             study_start_year=STUDY_START_YEAR_,
             study_end_year=STUDY_END_YEAR_,
             lifespan=AGE_MAX_IN_YEARS_,
             prg_intrvl=AGE_MOM_IN_YEARS_,
             bw_frac_lip=BODYWEIGHT_AND_LIPID_FRACTION_,
             timesteps_per_month=TIMESTEPS_PER_MONTH,
             average_lact_time=AVERAGE_LACT_TIME_,
             k_lac=K_LAC_GEN,
             k_elim=K_ELIM,
             peak_dose_intensity=PEAK_DOSING_INTENSITY_):
    t_start = study_start_year
    t_final = study_end_year
    delta_t = 1 / timesteps_per_month

    num_steps = np.floor((t_final - t_start) / delta_t) + 1

    if bm_data:
        bm_data_eval(bm_data,
                     cte,
                     kinetics_order,
                     plot_kinetics)

    cong_spag_df = pd.read_csv(cong_spag)
    congener2eval = bm_data_err(cte, cong_spag_df)

    intakeCall = asym_int_dist(peak_dose_intensity, 35, t_start, t_final, delta_t)
    # intakeCall = intake_int_dist(peak_dose_intensity, 70, t_start, t_final, delta_t)

    y = []

    (y, t) = generation_mass_balance(y=y,
                                     congener='hcb',
                                     gens=gens,
                                     simulation_start=t_start,
                                     simulation_end=t_final,
                                     delta_t=delta_t,
                                     intakeCall=intakeCall,
                                     prg_intrvl=prg_intrvl,
                                     lifespan=lifespan,
                                     num_steps=num_steps,
                                     bw_frac_lip=bw_frac_lip,
                                     k_lac=k_lac,
                                     average_lact_time=average_lact_time,
                                     k_elim=k_elim)


if __name__ == "__main__":
    pk_model(gens = 4)


