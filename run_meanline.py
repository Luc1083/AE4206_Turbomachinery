import meanline as ml

lieblein_model = ml.Lieblein_Model()
lieblein_model.plot_graphs(beta_start=0, beta_end=80, solidity_start=0.2, solidity_end=2.2, no_points_beta=20,
                           delta_solidity=0.2, tc_start=0, tc_end=0.2, no_points_tc=20)


fan = ml.Fan(Mach_inlet=0.6, AR_rotor=10, AR_stator=10, taper_rotor=0.7, taper_stator=0.5, n=1, no_blades_rotor=30,
             no_blades_stator=30, beta_tt=1.6, P0_cruise=39513.14, T0_cruise=250.13, mdot=80, omega=5000,
             hub_tip_ratio=0.3, gamma=1.4, R_air=287, eta_tt_estimated=0.9, Cp_air=1006, Cv_air=715.9,
             row_chord_spacing_ratio=0.5, lieblein_model=lieblein_model, profile="NACA-65")
ml.Fan_Plots(fan)