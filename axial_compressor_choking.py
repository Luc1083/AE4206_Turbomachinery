import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt


# Estimation of axial compressor choking point
# Author: ir. Andrea Giuffre'
# Reference: C. Freeman, N.A. Cumpsty - "A Method for the Prediction of Supersonic Compressor Blade Performance", 1989
# Delft University of Technology - All rights reserved


# compressor design data
beta_blade = np.array([])       # hub, mid, tip
R_hub = 0.3
R_tip = 0.6
t_le = 0.05
Nbl = 16

# operating conditions
Pt = 39513.14
Tt = 250.13
mass_flow_vec = np.arange(70, 90, 5)     # range of mass flow to analyze
rpm = 5000

# fluid information (perfect gas)
gamma = 1.4
R = 287.05
cp = 1004.5


def compute_inlet_flow_conditions(mass_flow, R_loc):
    """
    Given the operating point specified in terms of mass flow rate and rotational speed and
    the inlet geometry specified in terms of blade height, mean radius and blade angle,
    compute the inlet relative Mach number and the flow angle.
    """
    Area = 2 * np.pi * R_loc * H_in
    omega = rpm * (2 * np.pi) / 60
    U = omega * R_loc
    Mach_abs = opt.fsolve(massflow_isentropic, (0.5), args=(Area, mass_flow), full_output=False)
    T = Tt / (1 + (gamma - 1) / 2 * Mach_abs ** 2)
    P = Pt / ((Tt / T) ** (gamma / (gamma - 1)))
    Vm = Mach_abs * np.sqrt(gamma * R * T)
    Wm = Vm
    Wt = 0.0 - U
    W = np.sqrt(Wm ** 2 + Wt ** 2)
    Mach_rel = W / np.sqrt(gamma * R * T)
    beta = np.arctan(Wt / Wm)

    return Mach_rel, beta, T, P


def massflow_isentropic(p, *data):
    """
    Given the annulus area and the mass flow passing through it, compute the absolute Mach number, assuming
    isentropic flow and perfect gas.
    """
    A, m = data
    Mach = p

    m_computed = Pt * A / np.sqrt(R * Tt) * Mach * np.sqrt(gamma) * (1 + (gamma - 1) / 2 * Mach ** 2) ** \
                 ((1 + gamma) / (2 * (1 - gamma)))
    res = (m_computed - m) / m

    return res


def freeman_cv(t_th, beta, beta_blade, Mach_rel_in, Mach_rel_out):
    """
    Compute left-hand and right-hand sides of the equation derived by Freeman and Cumpsty for the estimation
    of choking point in axial compressor cascades, assuming perfect gas.
    Inlet and outlet refer to the boundaries of the control volume, not of the cascade.
    """
    lhs = ((1 + (gamma - 1) / 2 * (Mach_rel_out ** 2)) ** (- 1 / 2)) * \
          (1 + gamma * (Mach_rel_out ** 2) * (1 - t_th)) / (Mach_rel_out * (1 - t_th))
    rhs = ((1 + (gamma - 1) / 2 * (Mach_rel_in ** 2)) ** (- 1 / 2)) * \
          (np.cos(beta_blade) / np.cos(beta) + gamma * (Mach_rel_in ** 2) * np.cos(beta - beta_blade)) / Mach_rel_in
    
    return lhs, rhs


def compute_post_shocks_flow_conditions(p, *data):
    """
    Compute Mach number at the outlet of the control volume defined by Freeman and Cumpsty.
    """
    t_th, beta, beta_blade, Mach_rel_in = data
    Mach_rel_out = p

    lhs, rhs = freeman_cv(t_th, beta, beta_blade, Mach_rel_in, Mach_rel_out)
    res = (lhs - rhs) / lhs

    return res


def compute_entropy_generation(beta_blade, t_th, T_in, rho_in, P_in, Mach_rel_in, Mach_rel_out):
    """
    Given inlet state and outlet relative Mach number, compute entropy generation within the control volume,
    assuming perfect gas.
    """
    T_out = T_in * (1 + (gamma - 1) / 2 * Mach_rel_in ** 2) / (1 + (gamma - 1) / 2 * Mach_rel_out ** 2)
    rho_out = rho_in * (Mach_rel_in * np.sqrt(T_in) * np.cos(beta) / np.cos(beta_blade)) / \
              (Mach_rel_out * np.sqrt(T_out) * (1 - t_th))
    P_out = rho_out * R * T_out
    ds = cp * np.log(T_out / T_in) - R * np.log(P_out / P_in)

    return ds


# choose hub [0], mid [1] or tip [2]
loc = 1

# run
R_mean = (R_hub + R_tip) / 2
H_in = R_tip - R_hub
beta_blade = np.deg2rad(beta_blade)
Radius = np.array([R_hub, R_mean, R_tip])
throat = 2 * np.pi * Radius / Nbl * np.cos(beta_blade) - t_le
t_th = t_le / throat

incidence_vec = np.array([])
Mach_rel_in_vec = np.array([])
Mach_rel_out_vec = np.array([])
ds_vec = np.array([])
residual_vec = np.array([])

for mass_flow in mass_flow_vec:
    Mach_rel_in, beta, T_in, P_in = compute_inlet_flow_conditions(mass_flow, Radius[loc])
    incidence = np.rad2deg(np.abs(beta) - np.abs(beta_blade[loc]))
    result = opt.fsolve(compute_post_shocks_flow_conditions, (Mach_rel_in),
                        args=(t_th[loc], beta, beta_blade[loc], Mach_rel_in), full_output=True)
    Mach_rel_out = result[0]
    rho_in = P_in / (R * T_in)
    ds = compute_entropy_generation(beta_blade[loc], t_th[loc], T_in, rho_in, P_in, Mach_rel_in, Mach_rel_out)

    incidence_vec = np.append(incidence_vec, incidence)
    Mach_rel_in_vec = np.append(Mach_rel_in_vec, Mach_rel_in)
    Mach_rel_out_vec = np.append(Mach_rel_out_vec, Mach_rel_out)
    ds_vec = np.append(ds_vec, ds)
    residual_vec = np.append(residual_vec, result[1]['fvec'])

choking_idx = np.argmax(residual_vec > 1e-6) - 1

plt.figure()
plt.plot(mass_flow_vec, incidence_vec, 'black')
plt.plot(mass_flow_vec[choking_idx], incidence_vec[choking_idx], 'ro')
plt.xlabel('mass flow rate [kg/s]')
plt.ylabel('incidence angle [deg]')
plt.grid(1)

plt.figure()
plt.plot(mass_flow_vec, Mach_rel_in_vec, 'black')
plt.plot(mass_flow_vec[choking_idx], Mach_rel_in_vec[choking_idx], 'ro')
plt.xlabel('mass flow rate [kg/s]')
plt.ylabel('inlet Mach [-]')
plt.grid(1)

plt.figure()
plt.plot(mass_flow_vec, Mach_rel_out_vec, 'black')
plt.plot(mass_flow_vec[choking_idx], Mach_rel_out_vec[choking_idx], 'ro')
plt.xlabel('mass flow rate [kg/s]')
plt.ylabel('post shocks Mach [-]')
plt.grid(1)

plt.figure()
plt.plot(mass_flow_vec, ds_vec, 'black')
plt.plot(mass_flow_vec[choking_idx], ds_vec[choking_idx], 'ro')
plt.xlabel('mass flow rate [kg/s]')
plt.ylabel('entropy rise [J/(kg K)]')
plt.grid(1)

plt.figure()
plt.plot(mass_flow_vec, residual_vec, 'black')
plt.plot(mass_flow_vec[choking_idx], residual_vec[choking_idx], 'ro')
plt.xlabel('mass flow rate [kg/s]')
plt.ylabel('residual [-]')
plt.grid(1)
plt.show()
