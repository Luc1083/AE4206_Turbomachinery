import numpy as np


class FlowField:
    def __init__(self, ref_values, U, theta, r, Cp, gamma, rho, Vx, Vr, Vth, e):
        self.reference_values = ref_values
        self.U = U
        self.theta = theta
        self.r = r
        self.Cp = Cp
        self.gamma = gamma
        self.rho = rho
        self.Vx = Vx
        self.Vr = Vr
        self.Vth = Vth
        self.e = e

        self.variables = self.gen_variables()

    def gen_variables(self):
        var_dict = dict()

        U = self.U
        theta = self.theta
        Cp = self.Cp
        gamma = self.gamma
        Cv = Cp / gamma
        Rg = Cp - Cv
        rho = self.rho
        Vx = self.Vx
        Vr = self.Vr
        Vth = self.Vth
        e = self.e
        [N_i, N_j, N_k] = U.shape
        # var_dict['test'] = np.ones([N_i,N_j,N_k])

        # Absolute frame velocity
        var_dict['Vx_stn'] = Vx
        var_dict['Vy_stn'] = Vr * np.sin(theta) + Vth * np.cos(theta)
        var_dict['Vz_stn'] = Vr * np.cos(theta) - Vth * np.sin(theta)
        var_dict['Vr_stn'] = Vr
        var_dict['Vth_stn'] = Vth
        var_dict['Vm_stn'] = np.sqrt(Vx ** 2 + Vr ** 2)
        var_dict['pitch'] = np.arccos(Vx / var_dict['Vm_stn']) * np.sign(Vr) * 180 / np.pi
        var_dict['yaw_stn'] = np.arctan(Vth / var_dict['Vm_stn']) * 180 / np.pi
        var_dict['alpha'] = np.arctan(Vth / Vx) * 180 / np.pi

        # Relative frame velocity
        var_dict['U'] = U
        var_dict['Vx_rel'] = Vx
        var_dict['Vr_rel'] = Vr
        var_dict['Vth_rel'] = Vth - U
        var_dict['Vm_rel'] = var_dict['Vm_stn']
        var_dict['Vy_rel'] = Vr * np.sin(theta) + (Vth - U) * np.cos(theta)
        var_dict['Vz_rel'] = Vr * np.cos(theta) - (Vth - U) * np.sin(theta)
        var_dict['yaw_rel'] = np.arctan(var_dict['Vth_rel'] / var_dict['Vm_rel'])
        var_dict['beta'] = np.arctan(var_dict['Vth_rel'] / Vx) * 180 / np.pi

        # Thermodynamic properties
        var_dict['e_stn'] = e
        var_dict['k_stn'] = 0.5 * (Vx * Vx + Vr * Vr + Vth * Vth)
        var_dict['k_rel'] = 0.5 * (Vx * Vx + Vr * Vr + var_dict['Vth_rel'] * var_dict['Vth_rel'])
        var_dict['T'] = (e - var_dict['k_stn']) / Cv
        var_dict['rho'] = rho
        var_dict['P'] = Rg * rho * var_dict['T']

        var_dict['Tt_stn'] = var_dict['T'] + var_dict['k_stn'] / Cp
        var_dict['Pt_stn'] = var_dict['P'] * (var_dict['Tt_stn'] / var_dict['T']) ** (gamma / (gamma - 1))
        var_dict['rhot_stn'] = var_dict['Pt_stn'] / (Rg * var_dict['Tt_stn'])

        var_dict['Tt_rel'] = var_dict['T'] + var_dict['k_rel'] / Cp
        var_dict['Pt_rel'] = var_dict['P'] * (var_dict['Tt_rel'] / var_dict['T']) ** (gamma / (gamma - 1))
        var_dict['rhot_rel'] = var_dict['Pt_rel'] / (Rg * var_dict['Tt_rel'])

        Tt_ref = self.reference_values[0]
        Pt_ref = self.reference_values[1]
        var_dict['s'] = Cp * np.log(var_dict['T'] / Tt_ref) - Rg * np.log(var_dict['P'] / Pt_ref)

        # Mach number
        var_dict['a'] = np.sqrt(gamma * Rg * var_dict['T'])
        var_dict['V_stn'] = np.sqrt(2 * var_dict['k_stn'])
        var_dict['M_stn'] = var_dict['V_stn'] / var_dict['a']
        var_dict['V_rel'] = np.sqrt(2 * var_dict['k_rel'])
        var_dict['M_rel'] = var_dict['V_rel'] / var_dict['a']

        # # These are not really useful and they are easy enough to generate in ParaView: do not overload the .dat file
        # # Total to total efficiency
        var_dict['beta_tt'] = var_dict['Pt_stn'] / Pt_ref  # Total to total pressure ratio
        # Tt_is = Tt_ref*var_dict['beta_tt']**((gamma-1)/gamma) # Isentropic total temperature corresponding to the total pressure change
        # var_dict['w_is_tt'] = Cp*(Tt_is - Tt_ref)  # Isentropic work gained by the fluid (>0: compressor; <0: turbine)
        # var_dict['w_w_tt'] = Cp * (var_dict['Tt_stn'] - Tt_is) # Wasted work (always > 0, it is a work that the fluid always takes ('steals') from the machine)
        # var_dict['w_tt'] = Cp*(var_dict['Tt_stn'] - Tt_ref) # Fluid total work absorbed (total amount of energy absorbed by the flow, (>0: compressor; <0: turbine))
        #
        # # Total to static efficiency:
        var_dict['beta_ts'] = var_dict['P'] / Pt_ref  # Total to static pressure ratio
        # T_is = Tt_ref * var_dict['beta_ts'] ** ((gamma - 1) / gamma)  # Isentropic total temperature corresponding to the total to static pressure change
        # var_dict['w_is_ts'] = Cp * (T_is - Tt_ref)  # Isentropic work gained by the fluid (>0: compressor; <0: turbine)
        # var_dict['w_w_ts'] = Cp * (var_dict['T'] - T_is)  # Wasted work (always > 0, it is a work that the fluid always takes ('steals') from the machine)
        # var_dict['w_ts'] = var_dict['w_is_ts'] + var_dict['w_w_ts']  # Fluid total work absorbed (total amount of energy absorbed by the flow, (>0: compressor; <0: turbine))

        # Geometry useful to output as a variable:
        I = np.zeros([N_i, N_j, N_k])
        J = np.zeros([N_i, N_j, N_k])
        K = np.zeros([N_i, N_j, N_k])
        for i in range(N_i): I[i, :, :] = i
        for j in range(N_j): J[:, j, :] = j
        for k in range(N_k): K[:, :, k] = k

        var_dict['i_index'] = I / (N_i - 1)
        var_dict['j_index'] = J / (N_j - 1)
        var_dict['k_index'] = K / (N_k - 1)

        var_dict['r'] = self.r
        var_dict['theta'] = self.theta

        return var_dict

    # This was just a test, not really useful
    # def plot_contour_M(self):
    #     M = self.variables['T']
    #     M = M[:,:,36]
    #     plt.figure()
    #     plt.contourf(M)
    #     plt.xlabel('J')
    #     plt.ylabel('I')
    #     plt.colorbar()
