import meanline as ml
import numpy as np
import matplotlib.pyplot as plt
from pymoo.optimize import minimize
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.termination import get_termination

lieblein_model = ml.Lieblein_Model()
# lieblein_model.plot_graphs(beta_start=0, beta_end=80, solidity_start=0.2, solidity_end=2.2, no_points_beta=20,
                           # delta_solidity=0.2, tc_start=0, tc_end=0.2, no_points_tc=20)

# fan = ml.Fan(Mach_inlet=0.69996347, AR_rotor=2.93991714, AR_stator= 3.41924708, taper_rotor= 0.95974572, taper_stator= 0.71619419, n= 0.20165627, no_blades_rotor=  32,
#              no_blades_stator= 48, beta_tt=1.6, P0_cruise=39513.14, T0_cruise=250.13, mdot=80, omega=5000,
#              hub_tip_ratio= 0.39602997, gamma=1.4, R_air=287, eta_tt_estimated=0.9, row_chord_spacing_ratio=0.5, lieblein_model=lieblein_model,
#              profile="NACA-65", methodology="controlled vortex")
# #
#
#
# ml.Fan_Plots(fan)
#
# print('Theta: %s, Psi: %s, R: %s' % (fan.theta,fan.psi_mean, fan.R_mean))
# print('Stage Efficiency: ', fan.eta_tt_estimated)
# print('Mean Periferal Velocity U2: ', fan.U_mean)
# print('\nBeta 1 angle at leading edge of rotor: \n', np.degrees(fan.beta_1_rotor_distribution))
# print('\nSolidity Rotor Distribution: \n', fan.solidity_rotor_distribution)
# print('\nChamber Angle Distribution: \n', fan.chamber_angle_rotor)
# print('\nMach Rotor: \n',fan.Mach_rotor)
# print('\nMach Stator: \n', fan.Mach_stator)
# print('\nIncidence Angles Rotor: \n', fan.i_rotor)
# print('\nIncidence Angles Stator: \n', fan.i_stator)
# print('\nDeviation Angles Rotor: \n', fan.delta_rotor)
# print('\nDeviation Angles Stator: \n', fan.delta_stator)

# Fm_r, Ft_r = fan.rotor_force()
# print('\nStator Force ')
# print("\nAxial / Meridional Force: \n", Fm_r/1000)
# print("\nTangential Force: \n", Ft_r/1000)
#
# Fm_s, Ft_s = fan.stator_force()
# print('\nStator Force ')
# print("\nAxial / Meridional Force: \n", Fm_s/1000)
# print("\nTangential Force: \n", Ft_s/1000 )

#ml.Fan_Plots(fan)

fan_optimisation_problem = ml.optimize_design_elementwise()

algorithm = NSGA2(pop_size=50)
termination = get_termination("n_gen", 100)

result = minimize(problem=fan_optimisation_problem,
                  algorithm=algorithm,
                  termination=termination,
                  seed=1,
                  verbose=True)

X = result.X
F = result.F
G = result.G
print(X)
print(F)
print(G)

def write_array_to_file(file_path, array):
    try:
        with open(file_path, 'w') as file:
            for item in array:
                file.write(str(item) + '\n')
        print("Array successfully written to", file_path)
    except Exception as e:
        print("Error occurred while writing to file:", str(e))


file_path_1 = "Variable_Optimisation_Output.txt"
file_path_2 = "Objective_Optimisation_Output.txt"
write_array_to_file(file_path_1, X)
write_array_to_file(file_path_2, F)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(1-F[:, 0], F[:, 1], F[:, 2])
plt.show()
