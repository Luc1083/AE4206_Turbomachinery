import numpy as np
import scipy.optimize as opt


def optimize_design():
    ...


class Fan:
    def __init__(self, Mach_inlet, AR, taper, n, no_blades, beta_tt, P0_cruise, T0_cruise, PR, mdot, omega,
                 hub_tip_ratio, gamma, R_air):
        self.Mach_inlet = Mach_inlet
        self.AR = AR
        self.taper = taper
        self.n = n
        self.no_blades = no_blades
        self.betta_tt = beta_tt
        self.P0_cruise = P0_cruise
        self.T0_cruise = T0_cruise
        self.PR = PR
        self.mdot = mdot
        self.omega = omega

        self.gamma = gamma
        self.P_inlet = self.P0_cruise * (1 + Mach_inlet**2 * (gamma - 1) / 2)**(-gamma / (gamma - 1))
        self.T_inlet = self.T0_cruise * (1 + Mach_inlet**2 * (gamma - 1) / 2)**(-1)
        self.R_air = R_air
        self.rho_inlet = self.P_inlet / (self.R_air * self.T_inlet)
        self.Mach_inlet = Mach_inlet
        self.v_inlet = Mach_inlet * np.sqrt(gamma * R_air * self.T_inlet)

        self.A_ax = self.mdot / (self.v_inlet * self.rho_inlet)
        self.r_t_inlet = np.sqrt(self.A_ax / (np.pi * (1 - hub_tip_ratio**2)))
        self.r_t_hub = hub_tip_ratio * self.r_t_inlet
        self.r_mean = ...
        self.s_mean = ...
        self.h = ...
        self.c_mean = ...
        self.c_hub = ...
        self.c_tip = ...

        self.U = self.r_mean * omega * 2 * np.pi / 60
        self.theta = self.v_inlet / self.U

        # Calculate gamma, intermediate properties
        # Calculate angles, iterate
        # Calculate distribution

        # Calculate losses

        # Iterate until convergence


        # Plot function
        # Optimization function




    def estimate_efficiency(self):
        ...

    def generate_design(self):
        ...

    def establish_velocity_triangles(self):
        ...

    def calc_mixing_loss(self):
        ...

    def calc_BL_loss(self):
        ...

    def calc_endwall_loss(self):
        ...

    def calc_tip_leakage_loss(self):
        ...

    def calc_shock_loss(self):
        ...

    def calc_stall_margin(self):
        ...
