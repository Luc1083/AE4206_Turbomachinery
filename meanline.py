import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as intrp


def optimize_design():
    ...


class Fan:
    def __init__(self, Mach_inlet, AR_rotor, AR_stator, taper_rotor, taper_stator, n, no_blades_rotor, no_blades_stator,
                 beta_tt, P0_cruise, T0_cruise, mdot, omega, hub_tip_ratio, gamma, R_air, eta_tt_estimated, Cp_air,
                 Cv_air, row_chord_spacing_ratio, lieblein_model, profile):
        # Assign properties
        self.Mach_inlet = Mach_inlet
        self.AR_rotor = AR_rotor
        self.AR_stator = AR_stator
        self.taper_rotor = taper_rotor
        self.taper_stator = taper_stator
        self.n = n
        self.no_blades_rotor = no_blades_rotor
        self.no_blades_stator = no_blades_stator
        self.betta_tt = beta_tt
        self.P0_inlet = P0_cruise
        self.T0_inlet = T0_cruise
        self.mdot = mdot
        self.omega = omega
        self.hub_tip_ratio = hub_tip_ratio
        self.row_chord_spacing_ratio = row_chord_spacing_ratio
        self.lieblein_model = lieblein_model
        self.no_points = 50
        self.profile = profile

        # Calculate inlet properties
        self.gamma = gamma
        self.Cp = Cp_air
        self.Cv = Cv_air
        self.P_inlet = self.P0_inlet * (1 + (Mach_inlet**2) * ((gamma - 1) / 2)) ** (-gamma / (gamma - 1))
        self.T_inlet = self.T0_inlet * (1 + (Mach_inlet**2) * ((gamma - 1) / 2)) ** (-1)
        self.R_air = R_air
        self.rho_inlet = self.P_inlet / (self.R_air * self.T_inlet)
        self.Mach_inlet = Mach_inlet
        self.v_axial = Mach_inlet * np.sqrt(gamma * R_air * self.T_inlet)

        # Calculate inlet and rotor geometry
        self.A = self.mdot / (self.v_axial * self.rho_inlet)
        self.r_tip = np.sqrt(self.A / (np.pi * (1 - hub_tip_ratio ** 2)))
        self.r_hub_inlet_rotor = hub_tip_ratio * self.r_tip
        self.r_mean_rotor = np.sqrt(self.A / (2 * np.pi) + self.r_hub_inlet_rotor ** 2)
        self.h_blade_rotor = self.r_tip - self.r_hub_inlet_rotor
        self.s_mean_rotor = self.r_mean_rotor * 2 * np.pi / self.no_blades_rotor
        self.c_average_rotor = self.h_blade_rotor / self.AR_rotor
        self.c_hub_rotor = 2 * self.c_average_rotor / (1 + self.taper_rotor)
        self.c_tip_rotor = self.taper_rotor * self.c_hub_rotor
        self.c_mean_rotor = (((self.r_mean_rotor - self.r_hub_inlet_rotor) / self.h_blade_rotor) *
                             (self.c_tip_rotor - self.c_hub_rotor) + self.c_hub_rotor)
        self.rotor_solidity_mean = self.c_mean_rotor / self.s_mean_rotor

        # Calculate theta
        self.U_mean = self.r_mean_rotor * omega * 2 * np.pi / 60
        self.theta = self.v_axial / self.U_mean

        # Calculate intermediate total properties
        self.P0_exit_rotor = self.betta_tt * self.P0_inlet
        self.T0s_exit_rotor = self.T0_inlet * self.betta_tt ** ((self.gamma - 1) / self.gamma)
        self.T0_exit_rotor = (self.T0s_exit_rotor - self.T0_inlet) / eta_tt_estimated + self.T0_inlet
        self.psi_mean = (self.T0s_exit_rotor - self.T0_inlet) * self.Cp / (self.U_mean**2)

        # Calculate angles, velocities at the meanline
        self.beta_1 = np.arctan(- 1 / self.theta)
        self.R_mean = - self.psi_mean / 2 + 1
        self.beta_2 = np.arctan((self.psi_mean - 1) / self.theta)
        self.alpha_2 = np.arctan(np.tan(self.beta_2) + 1 / self.theta)
        self.w_1 = self.v_axial / np.cos(self.beta_1)
        self.v_2 = self.v_axial / np.cos(self.alpha_2)
        self.w_2 = self.v_axial / np.cos(self.beta_2)

        # Calculate intermediate static properties after rotor
        [self.T_exit_rotor, self.M_exit_rotor, self.P_exit_rotor, self.rho_exit_rotor] = (
            self.calc_intermediate_properties(self.v_2, self.T0_exit_rotor, self.P0_exit_rotor))

        # Create inlet geometry of stator
        self.A_exit_rotor = self.mdot / (self.v_axial * self.rho_exit_rotor)
        [self.r_hub_inlet_stator, self.r_mean_inlet_stator, self.h_blade_inlet_stator, self.s_mean_stator,
         self.c_mean_stator, self.c_hub_stator, self.c_tip_stator, self.stator_solidity_mean] = (
            self.calc_stator_geometry(self.A_exit_rotor))

        # Calculate intermediate static properties after stator
        [self.T_exit_stator, self.M_exit_stator, self.P_exit_stator, self.rho_exit_stator] = (
            self.calc_intermediate_properties(self.v_axial, self.T0_exit_rotor, self.P0_exit_rotor))

        # Create outlet geometry of stator
        self.A_exit_stator = self.mdot / (self.v_axial * self.rho_exit_stator)
        [self.r_hub_outlet_stator, self.r_mean_outlet_stator, self.h_blade_outlet_stator] \
            = self.calc_stator_geometry(self.A_exit_stator)[0:3]

        # Create all intermediate geometry, first for rotor
        self.x_rotor_hub_inlet = 0
        c_dummy = (((self.r_hub_inlet_stator - self.r_hub_inlet_rotor) / self.h_blade_rotor) *
                   (self.c_tip_rotor - self.c_hub_rotor) + self.c_hub_rotor)
        self.x_rotor_hub_outlet = self.c_hub_rotor - (self.c_hub_rotor - c_dummy) / 2
        self.x_rotor_tip_inlet = (self.c_hub_rotor - self.c_tip_rotor) / 2
        self.x_rotor_tip_outlet = self.c_hub_rotor - (self.c_hub_rotor - self.c_tip_rotor) / 2
        self.spacing_rows = c_dummy / self.row_chord_spacing_ratio

        # Then for stator
        self.x_stator_hub_inlet = self.x_rotor_hub_outlet + self.spacing_rows
        c_dummy = (((self.r_hub_outlet_stator - self.r_hub_inlet_stator) / self.h_blade_inlet_stator) *
                   (self.c_tip_stator - self.c_hub_stator) + self.c_hub_stator)
        self.x_stator_hub_outlet = self.c_hub_stator - (self.c_hub_stator - c_dummy) / 2 + self.x_stator_hub_inlet
        self.x_stator_tip_inlet = (self.c_hub_stator - self.c_tip_stator) / 2 + self.x_stator_hub_inlet
        self.x_stator_tip_outlet = (self.x_stator_hub_inlet + self.c_hub_stator -
                                    (self.c_hub_stator - self.c_tip_stator) / 2)

        # Find parameters for general whirl/ variable vortex for rotor
        self.vt_2 = np.abs(self.v_2 * np.sin(self.alpha_2))
        self.b_rotor = self.r_mean_rotor * self.vt_2 / 2
        self.a_rotor = self.vt_2 / (2 * self.r_mean_rotor ** n)

        # Find parameters for general whirl/ variable vortex for stator
        self.b_stator = - self.v_2 * np.sin(self.alpha_2) * self.r_mean_inlet_stator
        self.a_stator = self.v_2 * np.sin(self.alpha_2) / (2 * self.r_mean_inlet_stator ** 2)

        # Calculate distribution of angles and velocities for rotor
        self.r_rotor = np.linspace(self.r_hub_inlet_rotor, self.r_tip, self.no_points)
        [self.alpha_1_rotor_distribution, self.beta_1_rotor_distribution, self.U_rotor_distribution,
         self.w_1_rotor_distribution, self.v_2_rotor_distribution, self.alpha_2_rotor_distribution,
         self.w_2_rotor_distribution, self.beta_2_rotor_distribution, self.theta_rotor_distribution,
         self.psi_rotor_distribution, self.R_rotor_distribution, self.solidity_rotor_distribution,
         self.DF_rotor_distribution] = self.calc_angle_distribution_rotor(self.a_rotor, self.b_rotor, self.r_rotor)

        # Calculate distribution of angles and velocities for stator
        self.r_stator = np.linspace(self.r_hub_inlet_stator, self.r_tip, self.no_points)
        [self.alpha_2_stator_distribution, self.v_2_stator_distribution, self.alpha_1_stator_distribution,
         self.solidity_stator_distribution, self.DF_stator_distribution] = (
            self.calc_angle_distribution_stator(self.a_stator, self.b_stator, self.r_stator))

        # Calculate blade geometry, flow deviation and iterate

        # Calculate losses

        # Iterate until convergence

        # Calculate thickness

        # Plot function
        # Optimization function

    def calc_intermediate_properties(self, v_exit, T0_exit, P0_exit):
        v_exit_norm = v_exit / np.sqrt(self.Cp * T0_exit)
        M_exit = np.sqrt((1 / (self.gamma - 1)) * (v_exit_norm ** 2 / (1 - 0.5 * v_exit_norm ** 2)))
        P_exit = P0_exit * (1 + (M_exit ** 2) * ((self.gamma - 1) / 2)) ** (-self.gamma / (self.gamma - 1))
        T_exit = T0_exit * (1 + (M_exit ** 2) * ((self.gamma - 1) / 2)) ** (-1)
        rho_exit = P_exit / (self.R_air * T_exit)
        return [T_exit, M_exit, P_exit, rho_exit]

    def calc_stator_geometry(self, A):
        r_hub = np.sqrt(self.r_tip**2 - A / np.pi)
        r_mean = np.sqrt(A / (2 * np.pi) + r_hub**2)
        h_blade = self.r_tip - r_hub
        s_mean = r_mean * 2 * np.pi / self.no_blades_stator
        c_average = h_blade / self.AR_stator
        c_hub = 2 * c_average / (1 + self.taper_stator)
        c_tip = self.taper_stator * c_hub
        c_mean = ((r_mean - r_hub) / h_blade) * (c_tip - c_hub) + c_hub
        stator_solidity_mean = c_mean / s_mean
        return [r_hub, r_mean, h_blade, s_mean, c_mean, c_hub, c_tip, stator_solidity_mean]

    def calc_angle_distribution_rotor(self, a, b, r):
        alpha_1 = np.zeros(np.shape(r)[0])
        # Get U
        U = r * self.omega * 2 * np.pi / 60
        # Get w_1, beta_1
        w_1 = np.sqrt(self.v_axial ** 2 + U ** 2)
        beta_1 = np.arctan(-U / self.v_axial)
        # Get v_2, alpha_2
        vt_2 = a * r ** self.n + b / r
        v_2 = np.sqrt(vt_2 ** 2 + self.v_axial ** 2)
        alpha_2 = np.arctan(vt_2 / self.v_axial)
        # Get w_2, beta_2
        wt_2 = vt_2 - U
        w_2 = np.sqrt(wt_2 ** 2 + self.v_axial ** 2)
        beta_2 = np.arctan(wt_2 / self.v_axial)
        # Get R, psi, theta
        theta = 1 / (-np.tan(beta_1))
        psi = theta * np.tan(beta_2) + 1
        R = - psi / 2 + 1
        # Calculate DF and solidity
        wt_1 = np.abs(w_1 * np.sin(beta_1))
        solidity = np.linspace(self.c_hub_rotor, self.c_tip_rotor, self.no_points) * self.no_blades_rotor / (2 * np.pi * r)
        DF = (w_1 - w_2) / w_1 + (wt_1 - wt_2) / (2 * solidity * w_1)
        return [alpha_1, beta_1, U, w_1, v_2, alpha_2, w_2, beta_2, theta, psi, R, solidity, DF]

    def calc_angle_distribution_stator(self, a, b, r):
        # Calculate angles and velocities
        alpha_2 = np.zeros(np.shape(r)[0])
        vt_inlet = a * r ** self.n - b / r
        v_inlet = np.sqrt(self.v_axial ** 2 + vt_inlet ** 2)
        alpha_1 = np.arctan(vt_inlet / self.v_axial)
        # Calculate DF
        solidity = np.linspace(self.c_hub_stator, self.c_tip_stator, self.no_points) * self.no_blades_stator / (2 * np.pi * r)
        DF = (v_inlet - self.v_axial) / v_inlet + vt_inlet / (2 * solidity * v_inlet)
        return [alpha_2, v_inlet, alpha_1, solidity, DF]

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


class Fan_Plots:
    def __init__(self, Fan):
        self.Fan = Fan
        self.plot_thermodynamic_properties()
        self.plot_meridional_shape()
        self.plot_distribution_at_blade()
        self.plot_velocity_triangle(beta_1=self.Fan.beta_1, beta_2=self.Fan.beta_2, alpha_2=self.Fan.alpha_2,
                                    U=self.Fan.U_mean, location="meanline")
        self.plot_velocity_triangle(beta_1=self.Fan.beta_1_rotor_distribution[0],
                                    beta_2=self.Fan.beta_2_rotor_distribution[0],
                                    alpha_2=self.Fan.alpha_2_rotor_distribution[0],
                                    U=self.Fan.U_rotor_distribution[0], location="hub")
        self.plot_velocity_triangle(beta_1=self.Fan.beta_1_rotor_distribution[-1],
                                    beta_2=self.Fan.beta_2_rotor_distribution[-1],
                                    alpha_2=self.Fan.alpha_2_rotor_distribution[-1],
                                    U=self.Fan.U_rotor_distribution[-1], location="tip")
        self.print_data()

    def plot_thermodynamic_properties(self):
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()

        # Plot temperatures
        ax1.plot([self.Fan.x_rotor_hub_inlet, (self.Fan.x_rotor_hub_outlet + self.Fan.x_stator_hub_inlet) / 2,
                  self.Fan.x_stator_hub_outlet], [self.Fan.T0_inlet, self.Fan.T0_exit_rotor, self.Fan.T0_exit_rotor], "r",
                 label="Total temperature")
        ax1.plot([self.Fan.x_rotor_hub_inlet, (self.Fan.x_rotor_hub_outlet + self.Fan.x_stator_hub_inlet) / 2,
                  self.Fan.x_stator_hub_outlet], [self.Fan.T_inlet, self.Fan.T_exit_rotor, self.Fan.T_exit_stator], "c",
                 label="Static temperature")
        ax1.set_xlabel("Axial length (m)")
        ax1.set_ylabel("Temperature (K)")

        # Plot pressures
        ax2.plot([self.Fan.x_rotor_hub_inlet, (self.Fan.x_rotor_hub_outlet + self.Fan.x_stator_hub_inlet) / 2,
                  self.Fan.x_stator_hub_outlet], [self.Fan.P0_inlet, self.Fan.P0_exit_rotor, self.Fan.P0_exit_rotor],
                 label="Total pressure")
        ax2.plot([self.Fan.x_rotor_hub_inlet, (self.Fan.x_rotor_hub_outlet + self.Fan.x_stator_hub_inlet) / 2,
                  self.Fan.x_stator_hub_outlet], [self.Fan.P_inlet, self.Fan.P_exit_rotor, self.Fan.P_exit_stator],
                 label="Static pressure")
        ax2.set_ylabel("Pressure (Pa)")

        ax1.minorticks_on()
        ax1.grid(which='major', color='#DDDDDD', linewidth=0.8)
        ax1.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.8)
        ax1.legend(loc="upper left")
        ax2.legend(loc="lower right")
        ax2.minorticks_on()
        plt.savefig(f'thermodynamic_properties.png')
        plt.show()

    def plot_meridional_shape(self):
        # Plot rotor
        plt.plot([self.Fan.x_rotor_hub_inlet, self.Fan.x_rotor_hub_outlet, self.Fan.x_rotor_tip_outlet,
                  self.Fan.x_rotor_tip_inlet, self.Fan.x_rotor_hub_inlet],
                 [self.Fan.r_hub_inlet_rotor, self.Fan.r_hub_inlet_stator, self.Fan.r_tip,
                  self.Fan.r_tip, self.Fan.r_hub_inlet_rotor], label="Rotor")

        # Plot stator
        plt.plot([self.Fan.x_stator_hub_inlet, self.Fan.x_stator_hub_outlet, self.Fan.x_stator_tip_outlet,
                  self.Fan.x_stator_tip_inlet, self.Fan.x_stator_hub_inlet],
                 [self.Fan.r_hub_inlet_stator, self.Fan.r_hub_outlet_stator, self.Fan.r_tip, self.Fan.r_tip,
                  self.Fan.r_hub_inlet_stator], label="Stator")
        plt.xlabel("Axial length (m)")
        plt.ylabel("Radius (m)")
        plt.minorticks_on()
        plt.grid(which='major', color='#DDDDDD', linewidth=0.8)
        plt.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.8)
        plt.legend()
        plt.axis('scaled')
        plt.savefig(f'meridional_shape.png')
        plt.show()

    def plot_distribution_at_blade(self):
        fig = plt.figure()

        angles_graph = fig.add_subplot(211)
        angles_graph.plot(self.Fan.r_rotor, self.Fan.alpha_1_rotor_distribution * 180 / np.pi,
                          label="Rotor - alpha_1")
        angles_graph.plot(self.Fan.r_rotor, self.Fan.alpha_2_rotor_distribution * 180 / np.pi,
                          label="Rotor - alpha_2")
        angles_graph.plot(self.Fan.r_rotor, self.Fan.beta_1_rotor_distribution * 180 / np.pi,
                          label="Rotor - beta_1")
        angles_graph.plot(self.Fan.r_rotor, self.Fan.beta_2_rotor_distribution * 180 / np.pi,
                          label="Rotor - beta_2")
        angles_graph.plot(self.Fan.r_stator, self.Fan.alpha_1_stator_distribution * 180 / np.pi,
                          label="Stator - alpha_1")
        angles_graph.plot(self.Fan.r_stator, self.Fan.alpha_2_stator_distribution * 180 / np.pi,
                          label="Stator - alpha_2")
        angles_graph.set_xlabel("Radius (m)")
        angles_graph.set_ylabel("Angles (deg)")
        angles_graph.set_title("Distribution of angles along blade radius")

        angles_graph.minorticks_on()
        angles_graph.grid(which='major', color='#DDDDDD', linewidth=0.8)
        angles_graph.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.8)
        angles_graph.legend()

        properties_rotor_graph = fig.add_subplot(212)
        properties_rotor_graph.plot(self.Fan.r_rotor, self.Fan.theta_rotor_distribution, label="theta")
        properties_rotor_graph.plot(self.Fan.r_rotor, self.Fan.psi_rotor_distribution, label="psi")
        properties_rotor_graph.plot(self.Fan.r_rotor, self.Fan.R_rotor_distribution, label="R")
        properties_rotor_graph.plot(self.Fan.r_rotor, self.Fan.DF_rotor_distribution, label="DF - rotor")
        properties_rotor_graph.plot(self.Fan.r_stator, self.Fan.DF_stator_distribution, label="DF - stator")
        properties_rotor_graph.set_xlabel("Radius (m)")
        properties_rotor_graph.set_ylabel("Value (-)")
        properties_rotor_graph.set_title("Distribution of stage coefficients along rotor blade radius")

        properties_rotor_graph.minorticks_on()
        properties_rotor_graph.grid(which='major', color='#DDDDDD', linewidth=0.8)
        properties_rotor_graph.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.8)
        properties_rotor_graph.legend()

        plt.savefig(f'distribution_along_blade.png')
        plt.show()

    def plot_velocity_triangle(self, beta_1, beta_2, alpha_2, U, location):
        v_m = self.Fan.v_axial
        w_1_y = v_m * np.tan(beta_1)
        v_2_y = v_m * np.tan(alpha_2)
        w_2_y = v_m * np.tan(beta_2)

        # Plot the first triangle
        plt.arrow(0, 0, dx=v_m, dy=w_1_y, fc="blue",
                  head_width=6, width=2, length_includes_head=True, edgecolor="blue")
        plt.annotate(text="W_1", xy=(v_m / 2, w_1_y / 2))
        plt.arrow(0, 0, dx=v_m, dy=w_1_y + U, fc="green",
                  head_width=7, width=3, length_includes_head=True, edgecolor="green")
        plt.annotate(text="V_1", xy=(v_m / 2, (w_1_y + U) / 2))
        plt.arrow(v_m, w_1_y, dx=0, dy=U, fc="red",
                  head_width=4, head_length=16, width=1, length_includes_head=True, edgecolor="red")
        plt.annotate(text="U", xy=(v_m, w_1_y / 2))

        # Plot the second triangle
        plt.arrow(0, 0, dx=v_m, dy=v_2_y, fc="green",
                  head_width=6, width=2, length_includes_head=True, edgecolor="green")
        plt.annotate(text="V_2", xy=(v_m / 2, v_2_y / 2))
        plt.arrow(0, 0, dx=v_m, dy=w_2_y, fc="blue",
                  head_width=7, width=3, length_includes_head=True, edgecolor="blue")
        plt.annotate(text="W_2", xy=(v_m / 2, w_2_y / 2))
        plt.arrow(v_m, v_2_y - U, dx=0, dy=U, fc="red",
                  head_width=4, head_length=16, width=1, length_includes_head=True, edgecolor="red")
        plt.annotate(text="U", xy=(v_m, (v_2_y + w_2_y) / 2))

        plt.xlabel("V_x (m / s)")
        plt.ylabel("V_y (m / s)")
        plt.grid(which='major', color='#DDDDDD', linewidth=0.8)
        plt.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.8)
        plt.minorticks_on()
        plt.title(f"Velocity triangles at {location}")
        plt.savefig(f'velocity_triangle_{location}.png')
        plt.show()

    def print_data(self):
        print(f"Stage coefficients at meanline are: theta={self.Fan.theta}, psi={self.Fan.psi_mean}, "
              f"R={self.Fan.R_mean}")

        print(f"Angles at meanline are: alpha_1={0} deg, alpha_2={self.Fan.alpha_2 * 180 / np.pi} deg, "
              f"beta_1={self.Fan.beta_1 * 180 / np.pi} deg, beta_2={self.Fan.beta_2 * 180 / np.pi} deg")

        print(f"Velocities at meanline are: v1={self.Fan.v_axial} m/s, v2={self.Fan.v_2} m/s, w1={self.Fan.w_1} m/s, "
              f"w2={self.Fan.w_2} m/s ")

        print(f"Thermodynamic properties at meanline before rotor are: P0={self.Fan.P0_inlet} Pa, "
              f"T0={self.Fan.T0_inlet} K, P={self.Fan.P_inlet} Pa, T={self.Fan.T_inlet} K, M={self.Fan.Mach_inlet}")

        print(f"Thermodynamic properties at meanline after rotor are: P0={self.Fan.P0_exit_rotor}, "
              f"T0={self.Fan.T0_exit_rotor}, P={self.Fan.P_exit_rotor}, T={self.Fan.T_exit_rotor}, "
              f"M={self.Fan.M_exit_rotor}")

        print(f"Thermodynamic properties at meanline after stator are: P0={self.Fan.P0_exit_rotor} Pa, "
              f"T0={self.Fan.T0_exit_rotor} K, P={self.Fan.P_exit_stator} Pa, T={self.Fan.T_exit_stator} K, "
              f"M={self.Fan.M_exit_stator}")


class Lieblein_Model:
    def __init__(self):
        # Interpolate i_0 - beta_1 graph
        points, values = self.get_data_for_2Dinterpolation([("i0_solidity-0_4.csv", 0.4), ("i0_solidity-0_6.csv", 0.6),
                                                            ("i0_solidity-0_8.csv", 0.8), ("i0_solidity-1.csv", 1),
                                                            ("i0_solidity-1_2.csv", 1.2), ("i0_solidity-1_4.csv", 1.4),
                                                            ("i0_solidity-1_6.csv", 1.6), ("i0_solidity-1_8.csv", 1.8),
                                                            ("i0_solidity-2.csv", 2.0)], graph="i0")
        self.calc_i0 = intrp.LinearNDInterpolator(points=points, values=values)
        # Interpolate n
        points, values = self.get_data_for_2Dinterpolation([("n_solidity-0_4.csv", 0.4), ("n_solidity-0_6.csv", 0.6),
                                                            ("n_solidity-0_8.csv", 0.8), ("n_solidity-1.csv", 1),
                                                            ("n_solidity-1_2.csv", 1.2), ("n_solidity-1_4.csv", 1.4),
                                                            ("n_solidity-1_6.csv", 1.6), ("n_solidity-1_8.csv", 1.8),
                                                            ("n_solidity-2.csv", 2)], graph="n")
        self.calc_n = intrp.LinearNDInterpolator(points=points, values=values)
        # Interpolate Ki_t
        tc_array = np.genfromtxt(fname="lieblein_graphs/Ki_tc.csv", usecols=[0], delimiter=",", dtype=np.float32)
        Ki_t_array = np.genfromtxt(fname="lieblein_graphs/Ki_tc.csv", usecols=[1], delimiter=",", dtype=np.float32)
        self.calc_Ki_t = intrp.interp1d(x=tc_array, y=Ki_t_array, fill_value='extrapolate', kind="linear")

        # Interpolate delta
        points, values = self.get_data_for_2Dinterpolation([("sigma0_solidity-0_4.csv", 0.4),
                                                            ("sigma0_solidity-0_8.csv", 0.8),
                                                            ("sigma0_solidity-1.csv", 1),
                                                            ("sigma0_solidity-1.2.csv", 1.2),
                                                            ("sigma0_solidity-1_4.csv", 1.4),
                                                            ("sigma0_solidity-1_6.csv", 1.6),
                                                            ("sigma0_solidity-1_8.csv", 1.8),
                                                            ("sigma0_solidity-2.csv", 2)], graph="delta")
        self.calc_delta0 = intrp.LinearNDInterpolator(points=points, values=values)
        # Interpolate K_delta
        tc_array = np.genfromtxt(fname="lieblein_graphs/Ksigma_tc.csv", usecols=[0], delimiter=",", dtype=np.float32)
        K_delta_array = np.genfromtxt(fname="lieblein_graphs/Ksigma_tc.csv", usecols=[1], delimiter=",", dtype=np.float32)
        self.calc_K_delta = intrp.interp1d(x=tc_array, y=K_delta_array, fill_value='extrapolate', kind="linear")

        # Interpolate m coefficient for NACA-65
        beta1_array = np.genfromtxt(fname="lieblein_graphs/m_naca.csv", usecols=[0], delimiter=",", dtype=np.float32)
        m_array = np.genfromtxt(fname="lieblein_graphs/m_naca.csv", usecols=[1], delimiter=",", dtype=np.float32)
        self.calc_m_naca65 = intrp.interp1d(x=beta1_array, y=m_array, fill_value='extrapolate', kind="linear")

        # Interpolate m coefficient for DCA
        beta1_array = np.genfromtxt(fname="lieblein_graphs/m_dca.csv", usecols=[0], delimiter=",", dtype=np.float32)
        m_array = np.genfromtxt(fname="lieblein_graphs/m_dca.csv", usecols=[1], delimiter=",", dtype=np.float32)
        self.calc_m_dca = intrp.interp1d(x=beta1_array, y=m_array, fill_value='extrapolate', kind="linear")

        # Interpolate b coefficient
        beta1_array = np.genfromtxt(fname="lieblein_graphs/b_beta.csv", usecols=[0], delimiter=",", dtype=np.float32)
        b_array = np.genfromtxt(fname="lieblein_graphs/b_beta.csv", usecols=[1], delimiter=",", dtype=np.float32)
        self.calc_b = intrp.interp1d(x=beta1_array, y=b_array, fill_value='extrapolate', kind="linear")

    def get_data_for_2Dinterpolation(self, list, graph):
        beta_array = np.array([])
        solidity_array = np.array([])
        i0_array = np.array([])
        for [file, solidity] in list:
            file = "lieblein_graphs/" + file
            # Load data from files
            beta_dummy = np.genfromtxt(fname=file, usecols=[0], delimiter=",", dtype=np.float32)
            if not graph == "n":
                solidity_dummy = np.ones(np.shape(beta_dummy)[0] + 1) * solidity
            else:
                solidity_dummy = np.ones(np.shape(beta_dummy)[0]) * solidity
            i0_values = np.genfromtxt(fname=file, usecols=[1], delimiter=",", dtype=np.float32)
            # Append data to matrix
            if not graph == "n":
                beta_array = np.append(beta_array, [0])
                i0_array = np.append(i0_array, [0])
            beta_array = np.append(beta_array, beta_dummy)
            i0_array = np.append(i0_array, i0_values)
            solidity_array = np.append(solidity_array, solidity_dummy)
        points = np.transpose([beta_array, solidity_array])
        return points, i0_array

    def get_deviation_angle(self, beta_1, delta_beta, solidity, t_c, profile):
        # First calculate i_0
        if profile == "NACA-65":
            Ki_sh = 1.1
            Kdelta_sh = 1.1
            m = self.calc_m_naca65(beta_1)
        elif profile == "DCA":
            Ki_sh = 0.7
            Kdelta_sh = 0.75
            m = self.calc_m_dca(beta_1)
        Ki_t = self.calc_Ki_t(t_c)
        i0_10 = self.calc_i0(beta_1, solidity)
        i_0 = Ki_sh * i0_10 * Ki_t
        # Then calculate delta_0
        Kdelta_t = self.calc_K_delta(t_c)
        delta0_10 = self.calc_delta0(beta_1, solidity)
        delta_0 = Kdelta_sh * Kdelta_t * delta0_10
        # Calculate theta
        n = self.calc_n(beta_1, solidity)
        b = self.calc_b(beta_1)
        theta = (delta_beta + delta_0 - i_0) / (1 + n - m / (solidity ** b))  # camber line angle, not flow coefficient
        # Calculate optimal incidence and flow angle
        i = i_0 + n * theta
        delta = delta_0 + theta * m / (solidity ** b)
        return i, delta, theta

    def plot_graphs(self, beta_start, beta_end, solidity_start, solidity_end, no_points_beta, delta_solidity, tc_start,
                    tc_end, no_points_tc):
        # Create arrays of solidity and beta
        solidity_array = np.arange(solidity_start, solidity_end + delta_solidity, delta_solidity)
        beta_array = np.linspace(beta_start, beta_end, no_points_beta)
        tc_array = np.linspace(tc_start, tc_end, no_points_tc)
        # Create figures to store the data
        fig = plt.figure()

        i0_graph = fig.add_subplot(241)
        i0_graph.set_xlabel("Beta1 (deg)")
        i0_graph.set_ylabel("i0_10 (deg)")

        n_graph = fig.add_subplot(242)
        n_graph.set_xlabel("Beta1 (deg)")
        n_graph.set_ylabel("n coefficient (-)")

        Ki_t_graph = fig.add_subplot(243)
        Ki_t_graph.set_xlabel("max thickness (t/c)")
        Ki_t_graph.set_ylabel("coefficient Ki_t (-)")

        delta0_graph = fig.add_subplot(244)
        delta0_graph.set_xlabel("Beta1 (deg)")
        delta0_graph.set_ylabel("delta0_10 (deg)")

        Kdelta_t_graph = fig.add_subplot(245)
        Kdelta_t_graph.set_xlabel("max thickness (t/c)")
        Kdelta_t_graph.set_ylabel("Kdelta_t (-)")

        m_graph = fig.add_subplot(246)
        m_graph.set_xlabel("Beta1 (deg)")
        m_graph.set_ylabel("m coefficient (-)")

        b_graph = fig.add_subplot(247)
        b_graph.set_xlabel("Beta1 (deg)")
        b_graph.set_ylabel("Exponent b (-)")

        # Plot lines for constant solidity based on interpolation to verify it for all 2D graphs
        for solidity in solidity_array:
            dummy_solidity_array = np.ones(no_points_beta) * solidity
            # Plot i0
            i0_array = self.calc_i0(beta_array, dummy_solidity_array)
            i0_graph.plot(beta_array, i0_array, label=f"solidity={solidity:.2f}")
            # Plot n
            n_array = self.calc_n(beta_array, dummy_solidity_array)
            n_graph.plot(beta_array, n_array, label=f"solidity={solidity:.2f}")
            # Plot delta_0
            delta0_array = self.calc_delta0(beta_array, dummy_solidity_array)
            delta0_graph.plot(beta_array, delta0_array, label=f"solidity={solidity:.2f}")

        # Plot lines for 1D graphs
        Ki_t_graph.plot(tc_array, self.calc_Ki_t(tc_array))
        Kdelta_t_graph.plot(tc_array, self.calc_K_delta(tc_array))
        m_graph.plot(beta_array, self.calc_m_dca(beta_array), label="DCA")
        m_graph.plot(beta_array, self.calc_m_naca65(beta_array), label="NACA65")
        b_graph.plot(beta_array, self.calc_b(beta_array))

        for graph in [i0_graph, n_graph, Ki_t_graph, delta0_graph, Kdelta_t_graph, m_graph, b_graph]:
            graph.minorticks_on()
            graph.grid(which='major', color='#DDDDDD', linewidth=0.8)
            graph.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.8)
            # graph.legend()
        plt.savefig(f'interpolated_graphs.png')
        plt.show()






