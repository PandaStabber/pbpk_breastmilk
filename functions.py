# coding=utf-8
# !/usr/bin/env python3

import matplotlib.pyplot as plt
import scipy.io
import pandas as pd
from scipy import stats
import scipy.interpolate as interpolate
from scipy import integrate
from matplotlib.pylab import *


def generation_mass_balance(y,
                            congener,
                            gens,
                            simulation_start,
                            simulation_end,
                            delta_t,
                            intakeCall,
                            prg_intrvl,
                            lifespan,
                            num_steps,
                            bw_frac_lip,
                            k_lac,
                            average_lact_time,
                            k_elim):
    def body_mass(t, y):
        '''
        The first generation's mass balance should be specified above. Every other generations balance,
        assuming it's the same, can be specified here.

        Note that 'cntr' variable is a loop tracking tool. To specify that the previous box's mass balance
        should be employed, use itr_mtrx[0][cntr]-X), where X is the total number of mass balances - 1.
        This is because the array goes from 0 to 3, with length 4.

        Use the np.int() to surround the itr_mtrx calls because arrays should be tracked with
        integers, not floats.

        You can add more mass balances, but do not change dydt_matrix[0][gen] label.
        This is because these are a placeholder variables that are reorganized from an array to
        a matrix.

        For notes:
        aig_mother[gen]  # age the mother gives birth
        cbtg_child[gen]  # year child is born from previous gen
        aigd_mother[gen]  # age of death of the mother

       '''
        cntr = 0

        for gen in range(0, gens):

            k_lac_m2c_r = k_lac_m2c(t, gen, k_lac, cbtg_child, average_lact_time)
            k_lac_2cfm_r = k_lac_m2c(t, gen - 1, k_lac, cbtg_child, average_lact_time)

            if gen == 0:
                dydt_matrix[0][cntr] = bw_spl_der[0](t)
                dydt_matrix[1][cntr] = intakeCall(t) * y[0] \
                                       - k_elim * y[1] \
                                       - k_lac_m2c_r * y[1]
                cntr = np.int(cntr + 1)

            elif gen > 0:
                dydt_matrix[0][cntr] = bw_spl_der[gen](t)
                dydt_matrix[1][cntr] = intakeCall(t) * y[np.int(itr_mtrx[0][cntr])] \
                                       + k_lac_2cfm_r * y[np.int(itr_mtrx[1][cntr - 1])] \
                                       - k_lac_m2c_r * y[np.int(itr_mtrx[1][cntr])] \
                                       - k_elim * y[np.int(itr_mtrx[1][cntr])]

                cntr = np.int(cntr + 1)

        dydt = np.ravel(dydt_matrix, order='F')

        return dydt

    # dynamically set the number of odes as a function of gens
    n = len(y) * gens
    num_odes_in_gen = 2
    n = num_odes_in_gen * gens
    dydt = np.zeros((n, 1))

    # simulation bounds
    t_start = np.float(simulation_start)
    t_final = np.float(simulation_end)

    start_array_dydt = linspace(0, num_odes_in_gen * (gens - 1), gens)

    # automatically generate generational critical age definitions
    aig_mother = linspace(0,
                          prg_intrvl * (gens),
                          gens + 1)

    cbtg_child = linspace(prg_intrvl,
                          prg_intrvl * (gens + 1),
                          gens + 1)

    aigd_mother = linspace(lifespan,
                           (lifespan + (25 * gens)),
                           gens + 1)

    print("aig_mother",aig_mother)  # age the mother gives birth
    print("cbtg_child",cbtg_child)  # year child is born from previous gen
    print("aigd_mother",aigd_mother)  # age of death of the mother

    for gen in range(0, gens + 1):
        if np.all(gens >= 1):
            start_dydt_gen = []
            start_dydt_gen.append(start_array_dydt)
            for i in range(1, gens - 1):
                start_dydt_gen.append([x + i for x in start_array_dydt])
            start_dydt_gen = np.array(start_dydt_gen)

    odes_per_gen = range(0, num_odes_in_gen)
    dydt_matrix = np.zeros(shape=(len(odes_per_gen),
                                  gens),
                           dtype=object)

    order_array_counter = np.array(range(0, gens * len(odes_per_gen)))
    itr_mtrx = order_array_counter.reshape((len(odes_per_gen), gen),
                                           order='F')

    bw_spl_der = age_splines(gens,
                             aig_mother,
                             aigd_mother,
                             t_start,
                             t_final,
                             delta_t,
                             bw_frac_lip)

    t = np.zeros((np.int(num_steps), 1))

    # use ``vode`` with "backward differentiation formula" or 'bdf'
    r = integrate.ode(body_mass).set_integrator('vode',
                                                order=4,
                                                nsteps=num_steps,
                                                min_step=1e-10,
                                                method='bdf')

    y0 = np.zeros((np.int(gens * num_odes_in_gen), 1))
    r.set_initial_value(y0, t_start)

    # create vectors to store trajectories
    ode_init = np.zeros((np.int(num_steps) * gens * num_odes_in_gen))
    ode_init_matrix = ode_init.reshape((num_odes_in_gen * gens,
                                        np.int(num_steps)),
                                       order='F')

    iter_odes = range(0, num_odes_in_gen * gens, 1)

    # initialize k for while loop
    k = 1
    while r.successful() and k < num_steps:
        r.integrate(r.t + delta_t)
        t[k] = r.t
        for ode in iter_odes:
            ode_init_matrix[ode][k] = r.y[ode]
        k += 1

    for ode in iter_odes:
        ax1 = plt.subplot(len(iter_odes), 1, iter_odes[ode] + 1)
        plt.plot(t, ode_init_matrix[ode][:])
        ax1.plot(t, ode_init_matrix[ode][:])
        ax1.set_xlim(t_start, t_final)
        ax1.grid('on')

    plt.xlim(t_start, t_final)
    plt.legend(iter_odes)
    plt.show()

    return (y, t)


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
    return eval_bio_diff


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
                plt.plot(np.exp(x), np.exp(y_regression_line),
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
    # year_begin = 0
    num_steps_up = (peak_year - year_begin) / delta_t
    x_up_t_step = linspace(year_begin, peak_year, num=num_steps_up)

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

    # concatenate the up and down swing sides
    y_up_down = np.concatenate((upswing_reg_lin, dswing_reg_lin[1:]))
    x_up_down = np.concatenate((x_up_interp, x_d_interp[1:]))

    peak_intensity_spline = interpolate.InterpolatedUnivariateSpline(x_up_down, y_up_down)

    return peak_intensity_spline


def age_splines(gens, aig_mother, aigd_mother, year_begin, year_end, delta_t, bw_frac_lip):
    bw_frac_lip_df = pd.read_csv(bw_frac_lip)
    bw_frac_lip_df['mass_lipids_kg'] = bw_frac_lip_df.apply(cal_m_lip, 1)

    num_steps_before = []
    num_steps_after = []
    for gen in range(0, gens):
        num_steps_before.append(np.int((aig_mother[gen] - year_begin) / delta_t))
        num_steps_after.append(np.int((year_end - aigd_mother[gen]) / delta_t))

    age_spline_df = pd.DataFrame()
    for gen in range(0, gens):
        y_before = []
        x_before = []
        y_after = []
        x_after = []
        if num_steps_before[gen] == 0:
            y_gen = np.array(bw_frac_lip_df['mass_lipids_kg'])
            x_gen = linspace(aig_mother[gen], aigd_mother[gen], len(y_gen))

            if num_steps_after[gen] > 0:
                y_after = np.zeros(np.int(num_steps_after[gen]))
                x_after = linspace(aigd_mother[gen], year_end, num_steps_after[gen])

            # concatenate everything, but remove overlaps.
            y_all = np.concatenate([y_before[:-1], y_gen, y_after[1:-1]])
            x_all = np.concatenate([x_before[:-1], x_gen, x_after[1:-1]])

            age_spline = interpolate.InterpolatedUnivariateSpline(x_all, y_all).derivative()

        elif num_steps_before[gen] > 0:
            y_before = np.zeros(np.int(num_steps_before[gen]))
            x_before = linspace(aig_mother[gen - 1], aig_mother[gen], num_steps_before[gen])

            y_gen = np.array(bw_frac_lip_df['mass_lipids_kg'])
            x_gen = linspace(aig_mother[gen], aigd_mother[gen], len(y_gen))

            if num_steps_after[gen] > 0:
                y_after = np.zeros(np.int(num_steps_after[gen]))
                x_after = linspace(aigd_mother[gen], year_end, num_steps_after[gen])

            # concatenate everything, but remove overlaps.
            y_all = np.concatenate([y_before[:-1], y_gen, y_after[1:-1]])
            x_all = np.concatenate([x_before[:-1], x_gen, x_after[1:-1]])

            age_spline = interpolate.InterpolatedUnivariateSpline(x_all, y_all).derivative()

        age_spline_df = age_spline_df.append([age_spline], 1)

    return age_spline_df[0].ravel()

def k_lac_train(k_lac, gens, cbtg_child, aig_mother, aigd_mother, year_begin, year_end, delta_t, bw_frac_lip,
                average_lact_time):
    ''' look at the influence of variable lactation length. Variable lactation coefficients. '''
    num_steps_before = []
    num_steps_after = []
    gen_len = np.int(1e3)

    for gen in range(0, gens):
        num_steps_before.append(np.int((cbtg_child[gen] - year_begin) / delta_t))
        num_steps_after.append(np.int((year_end - (cbtg_child[gen] + average_lact_time[gen])) / delta_t))

    lac_spline_df = pd.DataFrame()

    for gen in range(0, gens):
        y_before = []
        x_before = []
        y_after = []
        x_after = []

        if num_steps_before[gen] == 0:
            y_gen = k_lac[gen] * gen_len
            x_gen = linspace(cbtg_child[gen], cbtg_child[gen] + average_lact_time[gen], gen_len)

            if num_steps_after[gen] > 0:
                y_after = np.zeros(np.int(num_steps_after[gen]))
                x_after = linspace(cbtg_child[gen] + average_lact_time[gen], year_end, num_steps_after[gen])

            # concatenate everything, but remove overlaps.
            y_all = np.concatenate([y_before[:-1], y_gen, y_after[1:-1]])
            x_all = np.concatenate([x_before[:-1], x_gen, x_after[1:-1]])

            lac_spline = interpolate.InterpolatedUnivariateSpline(x_all, y_all).derivative()

        elif num_steps_before[gen] > 0:
            y_before = np.zeros(np.int(num_steps_before[gen]))
            x_before = linspace(aig_mother[gen], cbtg_child[gen], num_steps_before[gen])

            y_gen = [k_lac[gen]] * gen_len
            x_gen = linspace(cbtg_child[gen], (cbtg_child[gen] + average_lact_time[gen]), gen_len)

            if num_steps_after[gen] > 0:
                y_after = np.zeros(np.int(num_steps_after[gen]))
                x_after = linspace((cbtg_child[gen] + average_lact_time[gen]), year_end, num_steps_after[gen])

            y_all = np.concatenate([y_before[:-1], y_gen, y_after[1:-1]])
            x_all = np.concatenate([x_before[:-1], x_gen, x_after[1:-1]])

            lac_spline = interpolate.UnivariateSpline(x_all, y_all)

        lac_spline_df = lac_spline_df.append([lac_spline], 1)

    return lac_spline_df[0].ravel()


def k_lac_m2c(t, gen, k_lac, cbtg_child, average_lact_time):
    # k_lac_r = 0.
    if (t >= cbtg_child[gen]) & (t <= cbtg_child[gen] + average_lact_time[gen]):
        k_lac_r = k_lac[gen]
    else:
        k_lac_r = 0.
    return k_lac_r


def bm_eval_(self, kin_order=1, plt_kin=True, out_kin=False, use_scipy_exp=False, p0=(1e-6, 1e-6, 1),
             use_scipy_lin2exp=False):
    kinetics_from_biomonitoring_data = []
    bm_data = pd.read_csv(self.bm_data)
    bm_data_t = bm_data.ix[:, 0]
    colors = plt.cm.rainbow(np.linspace(0, 1, len(self.cte)))

    def exp_fit(x, a, c, d):
        return a * np.exp(-c * x) + d

    x = np.array(bm_data_t)
    yr1 = x[0]

    for congener, c in zip(self.cte, colors):

        if kin_order == 1:

            x = x - x[0]
            y = np.array(bm_data[congener])

            if not all([use_scipy_exp, use_scipy_lin2exp]):

                if use_scipy_exp:

                    popt, pcov = curve_fit(exp_fit, x, y, p0=p0)

                    if out_kin:
                        print('congoner:', congener)
                        print('exp reg stats: ')
                        print("y = (", popt[0], ") * (e^(", -popt[1], " x )) + ", popt[2])
                        print('note: x can be a range of years (e.g., 2000-2010),')
                        print('or sequential (e.g.,''0-9). exp reg is performed by shifting x')
                        print('such that x[0] is zero. When plotted, it is shifted back.')

                    if plt_kin:
                        # shift the x values back to their respective years
                        x_shift = list(map(lambda x: x + yr1, x))
                        plt.plot(x_shift, y, color=c, marker='o', linestyle='',
                                 label=str('biomonitoring data; ' + congener))
                        plt.plot(x_shift, exp_fit(x, *popt),
                                 color=c, marker='', linestyle='-')

                if use_scipy_lin2exp:
                    y = np.log(np.array(bm_data[congener]))
                    (slope, intercept, r_value, p_value,
                     std_err) = stats.linregress(x, y)
                    y_regression_line = np.polyval([slope, intercept], x)

                    if out_kin:
                        print('congoner:', congener)
                        print('log-lin reg stats: ')
                        print("ln(y) = (", intercept, ") + (", slope, ") x")

                    if plt_kin:
                        # shift the x values back to their respective years
                        x_shift = list(map(lambda x: x + yr1, x))
                        plt.plot(x_shift, np.exp(y), color=c, marker='o', linestyle='',
                                 label=str('biomonitoring data; ' + congener))
                        plt.plot(x_shift, np.exp(y_regression_line),
                                 color=c, marker='', linestyle='-')


        elif kin_order == 2:
            x = np.log(np.array(bm_data_t))
            y = np.log(np.array(bm_data[congener]))
            (slope, intercept, r_value,
             p_value, std_err) = stats.linregress(x, y)

            y_regression_line = np.polyval([slope, intercept], x)
            kinetics_from_biomonitoring_data.append(slope)
            if plt_kin:
                plt.plot(np.exp(x), np.exp(y), color=c,
                         marker='o', linestyle='', label=str(
                        'biomonitoring data; ' + congener))
                plt.plot(np.exp(x), np.exp(y_regression_line),
                         color=c, marker='', linestyle='-')

    if plt_kin:
        plt.show()
