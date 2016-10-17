# coding=utf-8
import matplotlib.pyplot as plt
import scipy.io
import pandas as pd
from scipy import polyval, stats
import numpy as np
import os
import matplotlib.cm as cm
import scipy.interpolate as interpolate
from scipy.integrate import odeint
from scipy.ndimage.interpolation import shift
from scipy import integrate
from matplotlib.pylab import *


def mat_to_df(path):
    mat = scipy.io.loadmat(path,
                           matlab_compatible=True,
                           variable_names=None)
    mat_heads_init = list(mat)[0]
    mat_matrix = mat[mat_heads_init]
    df = pd.DataFrame(mat_matrix)
    return df


def bm_data_err(a, b):
    eval_bio_diff = list(set(a) - set(b))

    if not eval_bio_diff:
        print('''The following chemicals with biomonitoring data
         will not be evaluated'''.replace(
            '/t', ''), eval_bio_diff)
    else:
        print(''' Biomonitoring data may be available, but not
        evaluating chemicals:'''.replace(
            '/t', ''), eval_bio_diff)


# calculate mass of lipids for each month
def cal_m_lip(row):
    mass_lipids = row.avg_monthly_bodyweight_kg * row.lipid_fraction
    return mass_lipids


def bm_data_eval(bm_data, cte, kinetics_order, plot_kinetics):
    kinetics_from_biomonitoring_data = []
    bm_data = pd.read_csv(bm_data)
    bm_data_t = bm_data.ix[:, 0]
    bm_data_con = list(bm_data[-1:])
    bm_data_err(cte, bm_data_con)
    colors = cm.rainbow(np.linspace(0, 1, len(cte)))

    for congener, c in zip(cte, colors):

        if kinetics_order == 1:
            x = np.array(bm_data_t)
            y = np.log(np.array(bm_data[congener]))
            (slope, intercept, r_value, p_value,
             std_err) = stats.linregress(x, y)

            # TODO(output the kinetic rates and fit values to a table)
            y_regression_line = polyval([slope, intercept], x)
            kinetics_from_biomonitoring_data.append(slope)
            if plot_kinetics:
                plt.plot(x, np.exp(y), color=c, marker='o', linestyle='',
                         label=str('biomonitoring data; ' + congener))
                plt.plot(x, np.exp(y_regression_line),
                         color=c, marker='', linestyle='-')

        elif kinetics_order == 2:
            x = np.log(np.array(bm_data_t))
            y = np.log(np.array(bm_data[congener]))
            (slope, intercept, r_value,
             p_value, std_err) = stats.linregress(x, y)

            # TODO(output the kinetic rates and fit values to a table)
            y_regression_line = polyval([slope, intercept], x)
            kinetics_from_biomonitoring_data.append(slope)
            if plot_kinetics:
                plt.plot(np.exp(x), np.exp(y), color=c,
                         marker='o', linestyle='', label=str(
                        'biomonitoring data; ' + congener))
                plt.plot(np.exp(x), np.exp(yy_regression_line),
                         color=c, marker='', linestyle='-')

    if plot_kinetics:
        plt.show()


def intake_int_dist(peak_intensity, peak_year, year_begin, year_end, delta_t):
    # todo(remove peak_year,year_begin,year_end, and replace it with
    # default names in function)

    if (peak_year - year_begin) * 2 < (year_end - year_begin):
        print("unclipped symmetric intensity distribution")
        total_time_years = (year_end - year_begin)
        total_time_months = total_time_years * 12
        num_steps = total_time_months / delta_t

        # set total step space for simulation
        x_up_down = linspace(
            start=year_begin, stop=year_end, num=num_steps)

        upswing_time_years = peak_year - year_begin
        upswing_time_months = upswing_time_years * 12
        upswing_steps = upswing_time_months / delta_t
        # print upswing_steps

        downswing_steps = upswing_steps  # if symmetric

        # linearize the known datapoints and then create exponential
        # function:
        x_up = np.array([0, upswing_time_months])
        y_up = np.log([1e-10, peak_intensity])
        (slope_up, intercept_up, r_value, p_value,
         std_err) = stats.linregress(x_up, y_up)

        # interpolate using the provided slope and intercept to create the
        # regression line
        x_up_interp = linspace(0, upswing_time_months, upswing_steps)
        y_reg_line_up = polyval([slope_up, intercept_up], x_up_interp)
        # take the exp of the reg line to return it to an exp fit
        upswing_reg_lin = np.exp(y_reg_line_up)

        # downswing side of symmetric intensity distribution
        x_down = linspace(upswing_time_months,
                          upswing_time_months * 2, downswing_steps)
        y_down = upswing_reg_lin[::-1]  # reversed array

        # the remaining intensity is set to zero
        remaining_steps = num_steps - upswing_steps - downswing_steps
        y_remaining = np.zeros(np.int(remaining_steps), dtype=np.int)

        # concatenate the up and down swing sides
        y_up_down = np.concatenate((upswing_reg_lin, y_down, y_remaining))

        # create a callable spline function for the up and down functions
        peak_intensity_spline = interpolate.InterpolatedUnivariateSpline(
            x_up_down, y_up_down)

        return peak_intensity_spline

    else:
        print("clipped symettric intensity distribution")
        # todo(add symettric intensity distribution when clipped)


def asym_int_dist(peak_intensity, peak_year, year_begin, year_end, delta_t):
    year_begin = 0
    num_steps_up = (peak_year - year_begin) / delta_t
    x_up_t_step = linspace(year_begin, peak_year, num=num_steps_up)
    print(x_up_t_step)

    x_up = np.array([0, peak_year - year_begin])
    y_up = np.log([1e-10, peak_intensity])
    (slope_up, intercept_up, r_value, p_value,
     std_err) = stats.linregress(x_up, y_up)

    # interpolate using the provided slope and intercept to create the
    # regression line
    x_up_interp = linspace(0, peak_year, num_steps_up)
    y_reg_line_up = polyval([slope_up, intercept_up], x_up_interp)
    # take the exp of the reg line to return it to an exp fit
    upswing_reg_lin = np.exp(y_reg_line_up)

    # the remaining intensity is set to zero
    num_steps_d = (year_end - (peak_year + delta_t)) / delta_t
    x_d_interp = linspace(peak_year + delta_t, year_end, num_steps_d)
    x_d = np.array([peak_year + delta_t, year_end])
    y_d = np.log([peak_intensity, 1e-10])
    (slope_d, intercept_d, r_value, p_value, std_err) = stats.linregress(x_d, y_d)

    y_reg_line_down = polyval([slope_d, intercept_d], x_d_interp)
    dswing_reg_lin = np.exp(y_reg_line_down)
    print(dswing_reg_lin)

    # concatenate the up and down swing sides
    y_up_down = np.concatenate((upswing_reg_lin, dswing_reg_lin[1:]))
    x_up_down = np.concatenate((x_up_interp, x_d_interp[1:]))

    peak_intensity_spline = interpolate.InterpolatedUnivariateSpline(x_up_down, y_up_down)
    print(peak_intensity_spline(peak_year))

    return peak_intensity_spline


def age_splines(gens, aig_mother, aigd_mother, year_begin, year_end, delta_t, bw_frac_lip):
    bw_frac_lip_df = pd.read_csv(bw_frac_lip)
    bw_frac_lip_df['mass_lipids_kg'] = bw_frac_lip_df.apply(cal_m_lip, 1)

    num_steps_before = []
    num_steps_after = []
    for gen in range(0, gens):
        print(gen)
        num_steps_before.append(np.int((aig_mother[gen] - year_begin) / delta_t))
        num_steps_after.append(np.int((year_end - aigd_mother[gen]) / delta_t))

    age_spline_df = pd.DataFrame()
    for gen in range(0, gens):
        y_before = []
        x_before = []
        y_after = []
        x_after = []
        if num_steps_before[gen] == 0:

            # add in spline
            y_gen = np.array(bw_frac_lip_df['mass_lipids_kg'])
            x_gen = linspace(aig_mother[gen], aigd_mother[gen], len(y_gen))

            if num_steps_after[gen] > 0:
                # create zeros from death of gen 1 to the end
                y_after = np.zeros(np.int(num_steps_after[gen]))
                x_after = linspace(aigd_mother[gen], year_end, num_steps_after[gen])

            # concatenate everything, but remove overlaps.
            y_all = np.concatenate([y_before[:-1], y_gen, y_after[1:-1]])
            x_all = np.concatenate([x_before[:-1], x_gen, x_after[1:-1]])

            age_spline = interpolate.InterpolatedUnivariateSpline(x_all, y_all).derivative()

            age_spline_df = age_spline_df.append([age_spline], 1)


        elif num_steps_before[gen] > 0:
            y_before = np.zeros(np.int(num_steps_before[gen]))
            x_before = linspace(aig_mother[gen - 1], aig_mother[gen], num_steps_before[gen])

            # add in spline
            y_gen = np.array(bw_frac_lip_df['mass_lipids_kg'])
            x_gen = linspace(aig_mother[gen], aigd_mother[gen], len(y_gen))

            if num_steps_after[gen] > 0:
                # create zeros from death of gen 1 to the end
                y_after = np.zeros(np.int(num_steps_after[gen]))
                x_after = linspace(aigd_mother[gen], year_end, num_steps_after[gen])

            # concatenate everything, but remove overlaps.
            y_all = np.concatenate([y_before[:-1], y_gen, y_after[1:-1]])
            x_all = np.concatenate([x_before[:-1], x_gen, x_after[1:-1]])

            age_spline = interpolate.InterpolatedUnivariateSpline(x_all, y_all).derivative()

            age_spline_df = age_spline_df.append([age_spline], 1)

    return age_spline_df[0].ravel()


# todo(This function in progress)
def k_lac_train(gens, cbtg_child, aig_mother, aigd_mother, year_begin, year_end, delta_t):
    num_steps_before = []
    num_steps_after = []
    for gen in range(0, gens):
        print(gen)
        num_steps_before.append(np.int((cbtg_child[gen] - year_begin) / delta_t))
        num_steps_after.append(np.int((year_end - aigd_mother[gen]) / delta_t))

    print(num_steps_before)

    age_spline_df = pd.DataFrame()
    for gen in range(0, gens):
        y_before = []
        x_before = []
        y_after = []
        x_after = []
        if num_steps_before[gen] == 0:
            pass

            # add in spline
            # y_gen = np.array(bw_frac_lip_df['mass_lipids_kg'])
            # x_gen = linspace(aig_mother[gen], aigd_mother[gen], len(y_gen))
            #
            # if num_steps_after[gen] > 0:
            #     # create zeros from death of gen 1 to the end
            #     y_after = np.zeros(np.int(num_steps_after[gen]))
            #     x_after = linspace(aigd_mother[gen], year_end, num_steps_after[gen])
            #
            # # concatenate everything, but remove overlaps.
            # y_all = np.concatenate([y_before[:-1], y_gen, y_after[1:-1]])
            # x_all = np.concatenate([x_before[:-1], x_gen, x_after[1:-1]])
            #
            # age_spline = interpolate.InterpolatedUnivariateSpline(x_all, y_all).derivative()
            #
            # age_spline_df = age_spline_df.append([age_spline], 1)


            # elif num_steps_before[gen] > 0:
            #     y_before = np.zeros(np.int(num_steps_before[gen]))
            #     x_before = linspace(aig_mother[gen - 1], aig_mother[gen], num_steps_before[gen])

            # # add in spline
            # y_gen = np.array(bw_frac_lip_df['mass_lipids_kg'])
            # x_gen = linspace(aig_mother[gen], aigd_mother[gen], len(y_gen))
            #
            # if num_steps_after[gen] > 0:
            #     # create zeros from death of gen 1 to the end
            #     y_after = np.zeros(np.int(num_steps_after[gen]))
            #     x_after = linspace(aigd_mother[gen], year_end, num_steps_after[gen])
            #
            # # concatenate everything, but remove overlaps.
            # y_all = np.concatenate([y_before[:-1], y_gen, y_after[1:-1]])
            # x_all = np.concatenate([x_before[:-1], x_gen, x_after[1:-1]])
            #
            # age_spline = interpolate.InterpolatedUnivariateSpline(x_all, y_all).derivative()
            #
            # age_spline_df = age_spline_df.append([age_spline], 1)

            # return age_spline_df[0].ravel()


gens = 3
delta_t = 0.1
aig_mother = [0., 25., 50.]
year_begin = 0
aigd_mother = [80., 105., 130.]
year_end = 120
cbtg_child = [25., 50., 75.]

k_lac_train(gens, cbtg_child, aig_mother, aigd_mother, year_begin, year_end, delta_t)
