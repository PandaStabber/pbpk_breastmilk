def body_balance(t, y):

        # todo(put into function the k_elimination as a generic)
        n = len(y)  # 2: implies we have n ODEs
        dydt = np.zeros((n, 1))
        k_elimination = np.log(2) / (5 * 12) / 10

        def intake_amount_child(t, congener, mass_child):
            # todo(create automated peak year to run in loop based on the congener)
            year_begin = pop_start_year
            year_peak = 1977
            days_per_month = 365 / 12

            return_value = []
            if t < average_lact_time_months_:
                intake_curve_track_child.append(0)
                return 0

            elif t >= average_lact_time_months_ and t < (year_peak - year_begin) * 12:
                return_value = upswing_callable_intake(t) * mass_child * days_per_month / timesteps_per_month
                intake_curve_track_child.append(return_value)
                return upswing_callable_intake(t) * mass_child * days_per_month / timesteps_per_month

            elif t >= (year_peak - year_begin) * 12:  # and t < 2 * (year_peak - year_begin) * 12:

                return_value = upswing_callable_intake(t) * mass_child * days_per_month / timesteps_per_month
                intake_curve_track_child.append(return_value)
                return upswing_callable_intake(t) * mass_child * days_per_month / timesteps_per_month
            else:
                intake_curve_track_child.append(0)
                return 0

        def intake_amount(t, congener, mass_mother):
            # todo(create automated peak year to run in loop based on the congener)
            year_begin = pop_start_year
            year_peak = 1977
            days_per_month = 365 / 12

            # return_value = []
            if t < (year_peak - year_begin) * 12:
                return_value = upswing_callable_intake(t) * mass_mother * days_per_month / timesteps_per_month

            elif t >= year_peak * 12 and t < 2 * year_peak * 12:
                return_value = downswing_callable_intake(t) * mass_mother * days_per_month / timesteps_per_month

            else:
                return_value = 0

            intake_curve_track_mother.append(return_value)
            return return_value

        # mother mass, mass balance
        k_lactation = 0.005

        if generations == 1:
            if t < age_mom_in_years * 12:
                dydt[0] = spl_int_bw_lip_mass_deriv(t)
                dydt[1] = intake_amount(t, 'hcb', y[0]) - k_elimination * y[0]  # no lactation
                dydt[2] = 0
                dydt[3] = 0

            elif t >= age_mom_in_years * 12 and t < (age_mom_in_years * 12 + average_lact_time_months_):
                dydt[0] = spl_int_bw_lip_mass_deriv(t)
                dydt[1] = intake_amount(t, 'hcb', y[0]) - k_elimination * y[0] - k_lactation * y[0]
                dydt[2] = intake_amount_child(t, 'hcb', y[3]) - k_elimination * y[3] + k_lactation * y[0]
                dydt[3] = spl_int_bw_lip_mass_deriv(t - 12 * age_mom_in_years)

            elif t >= (age_mom_in_years * 12 + average_lact_time_months_) and t < (age_max_in_years * 12):
                dydt[0] = spl_int_bw_lip_mass_deriv(t)
                dydt[1] = intake_amount(t, 'hcb', y[0]) - k_elimination * y[0]
                dydt[2] = intake_amount_child(t, 'hcb', y[3]) - k_elimination * y[3] + k_lactation * y[0]
                dydt[3] = spl_int_bw_lip_mass_deriv(t - 12 * age_mom_in_years)

            elif t >= (age_max_in_years * 12):
                y[0] = 0

        if generations > 1:
            # todo(create function defining the number of generations, and the outputs for each y)
            print ('not quite ready for multi-generation use')

        # age_of_death_coefficient(t, y)
        return dydt

    # print intake_curve_track_child

    # use ``vode`` with "backward differentiation formula" or 'bdf'
    r = integrate.ode(body_balance).set_integrator('vode', order=5, nsteps=num_steps, max_step=delta_t, min_step=0,
                                                   method='bdf')
    r.set_initial_value([d_M_m_i, d_M_chem_m_i, d_C_chem_m_i, d_C_m_m_i], t_start)

    # create vectors to store trajectories
    t = np.zeros((np.int(num_steps), 1))

    mother_mass = np.zeros((np.int(num_steps), 1))
    mass_chemical_in_mother = np.zeros((np.int(num_steps), 1))
    mass_chemical_in_child = np.zeros((np.int(num_steps), 1))
    child_mass = np.zeros((np.int(num_steps), 1))