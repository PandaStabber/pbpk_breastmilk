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
        print '''The following chemicals with biomonitoring data
         will not be evaluated'''.replace(
            '/t', ''), eval_bio_diff
    else:
        print ''' Biomonitoring data may be available, but not
        evaluating chemicals:'''.replace(
            '/t', ''), eval_bio_diff


# calculate mass of lipids for each month
def cal_m_lip(row):
    mass_lipids = row.avg_monthly_bodyweight_kg * row.lipid_fraction
    return mass_lipids


def bm_data_eval(bm_data, cte, kinetics_order, plot_kinetics):

    kinetics_from_biomonitoring_data = []
    bm_data = pd.read_csv(bm_data)
    bm_data_t = bm_data.ix[:, 0]
    bm_data_con = list(bm_data[-1:])
    bm_data_err(
        cte, bm_data_con)
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
        y_up = np.log([1, peak_intensity])
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
