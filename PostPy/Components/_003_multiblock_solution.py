import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
from Components._004_grid import Coordinates
from Components._005_flow_field import FlowField


class BlockCFD:
    def __init__(self, data_grid, data_flow, reference_values):
        self.raw_data_grid = data_grid
        self.raw_data_flow = data_flow

        self.N_i = data_grid[1][0]
        self.N_j = data_grid[1][1]
        self.N_k = data_grid[1][2]

        self.grid = self.compute_grid()
        self.flow_field = self.compute_flow(reference_values)

    def compute_grid(self):
        # Size of the grid:
        N_i = self.N_i
        N_j = self.N_j
        N_k = self.N_k

        # Coordinates of the grid
        coordinates = self.raw_data_grid[6:]

        r = np.zeros([N_i, N_j, N_k])  # radial coordinate
        rth = np.zeros([N_i, N_j, N_k])  # tangential coordinate
        ax = np.zeros([N_i, N_j, N_k])  # meridional coordinate

        for j in range(N_j):
            for k in range(N_k):
                for i in range(N_i):
                    ax[i, j, k] = coordinates[k + j * N_k][0]
                    r[i, j, k] = coordinates[k + j * N_k][1]
                    rth[i, j, k] = coordinates[k + j * N_k][2 + i]

        return Coordinates(ax, r, rth)

    def compute_flow(self, ref_values):
        # Size of the grid:
        N_i = self.N_i
        N_j = self.N_j
        N_k = self.N_k

        # Extract standard variables:
        r = self.grid.r
        x = self.grid.x
        theta = self.grid.rth / r
        U = self.raw_data_grid[4][0] * self.grid.r
        Cp = self.raw_data_grid[2][0]
        gamma = self.raw_data_grid[2][1]
        rho = np.zeros([N_i, N_j, N_k])
        Vx = np.zeros([N_i, N_j, N_k])
        Vr = np.zeros([N_i, N_j, N_k])
        Vth = np.zeros([N_i, N_j, N_k])
        e = np.zeros([N_i, N_j, N_k])

        for k in range(N_k):
            for j in range(N_j):
                rho[:, j, k] = self.raw_data_flow[1][0 + j * N_i + k * (N_i * N_j):N_i + j * N_i + k * (N_i * N_j)]
                Vx[:, j, k] = self.raw_data_flow[2][0 + j * N_i + k * (N_i * N_j):N_i + j * N_i + k * (N_i * N_j)]
                Vr[:, j, k] = self.raw_data_flow[3][0 + j * N_i + k * (N_i * N_j):N_i + j * N_i + k * (N_i * N_j)]
                Vth[:, j, k] = self.raw_data_flow[4][0 + j * N_i + k * (N_i * N_j):N_i + j * N_i + k * (N_i * N_j)]
                e[:, j, k] = self.raw_data_flow[5][0 + j * N_i + k * (N_i * N_j):N_i + j * N_i + k * (N_i * N_j)]

        Vx = Vx / rho
        Vr = Vr / rho
        Vth = Vth / rho
        e = e / rho

        return FlowField(ref_values, U, theta, r, Cp, gamma, rho, Vx, Vr, Vth, e)

    def plot_Cp(self,k):
        N_i = self.N_i
        N_k = self.N_k

        k = int(np.float(k*(N_k-1)))
        P1 = self.flow_field.variables['P'][0,:,k]
        P2 = self.flow_field.variables['P'][N_i-1, :, k]
        Pt = self.flow_field.variables['Pt_rel'][0, 0, k]
        x = self.grid.x[0,:,k]

        Cp1 = (P1 - P1[0]) / (Pt - P1[0])
        Cp2 = (P2 - P1[0]) / (Pt - P1[0])

        plt.figure()
        plt.plot(x,Cp1)
        plt.plot(x,Cp2)
        plt.title(r'$C_p$ at grid surface K=' + str(k+1))
        plt.xlabel('x [m]')
        plt.ylabel(r'$C_p=\frac{(P_{x}-P_{1})}{(P_{1t,\ rel}-P_{1})}$')
        plt.show()

    def plot_variable_on_contour(self,variable,k):
        N_i = self.N_i
        N_k = self.N_k

        k = int(np.float(k*(N_k-1)))
        var1 = self.flow_field.variables[variable][0,:,k]
        var2 = self.flow_field.variables[variable][N_i-1, :, k]
        x = self.grid.x[0,:,k]

        plt.figure()
        plt.plot(x,var1)
        plt.plot(x, var2)
        plt.title(variable + ' at grid surface K=' + str(k+1))
        plt.xlabel('x [m]')
        plt.ylabel(variable + '[SI Units]]')
        plt.show()

    def plot_pitch_average_evolution(self, variable, avr_type='massFlowAve'):
        N_j = self.N_j
        plt.figure()
        for j in range(N_j):
            colour = (j / (N_j - 1), (N_j - 1 - j) / (N_j - 1), 0)
            self.plot_pitch_average_process(variable, j / (N_j - 1), avr_type, colour)
        plt.suptitle('Green: J=1;  Red: J=' + str(N_j))
        plt.title('Pitch-wise ' + avr_type + ' of ' + variable)

    def plot_pitch_average(self, variable, j, avr_type='massFlowAve'):
        plt.figure()
        self.plot_pitch_average_process(variable, j, avr_type)

    def plot_pitch_average_process(self, variable, j, avr_type='massFlowAve', colour=(0, 0, 0)):
        N_j = self.N_j
        J = int(np.round(j * (N_j - 1)))
        r = self.grid.r[0, J, :].tolist()
        avr = self.get_pitch_average(variable, j, avr_type).tolist()

        plt.plot(avr, r, color=colour)
        plt.xlabel(variable + ' [SI units]')
        plt.ylabel('r [m]')
        plt.title('Pitch-wise ' + avr_type + ' of ' + variable + ' at j=' + str(J + 1) + '/' + str(N_j))
        plt.show()

    def plot_B2B(self, variable, k, levels=None):
        plt.figure()
        self.plot_stream_surface_process(variable,k,levels)

    def plot_B2B_process(self, variable, k, levels=None):
        N_k = self.N_k
        k = int(np.round(k*(N_k-1)))

        var = self.flow_field.variables[variable][:,:,k]
        x = self.grid.x[:,:,k]
        rth = self.grid.rth[:,:,k]
        if levels == None:
            plt.contourf(x, rth, var)
        else:
            plt.contourf(x, rth, var,levels)
        plt.show()

    def get_pitch_average(self, variable: str, j: float, avr_type='massFlowAve') -> np.ndarray:
        """ Method to compute the pith-wise average of the fluid property 'variable' in an axial grid station J=j·N_j. It
            returns the span-wise distribution as a nparray. The method can perform different types of averages, and
            they are controlled by the parameter 'avr_type', which is set to mass flow average by default.

            Parameters
            ----------
            variable : str
                Variable that is going to be averaged. For mass flow and area averages it must be a key in the dictionary
                BlockCDF.flow_field.variables. If the average type is 'mixedOutAve' it can only take the values: 'P', 'Pt',
                'T', 'Tt', 'u_ax', 'u_th' or 's'. All of them in the stationary frame of reference.
            j : float
                Axial position where the average is performed as a fraction of the total domain length.
            avr_type : str
                Type of average to perform. it can be: 'massFlowAve', 'areaAve' or 'mixedOutAve'

            Returns
            -------
            avr : np.ndarray
               Distribution of the averaged quantity over the span-wise grid nodes.
        """

        avr = []
        N_k = self.N_k
        k_list = [[0.0], np.linspace(0, 1, N_k).tolist(), [1.0]]
        k_list = [item for sublist in k_list for item in sublist]

        if avr_type == 'areaAve':
            for i in range(1, N_k + 1):
                avr.append(self.get_area_average(variable, j, [0, 1], [k_list[i - 1], k_list[i + 1]]))
        elif avr_type == 'massFlowAve':
            for i in range(1, N_k + 1):
                avr.append(self.get_mass_flow_average(variable, j, [0, 1], [k_list[i - 1], k_list[i + 1]]))
        elif avr_type == 'mixedOutAve':
            variable_names = ['P', 'Pt', 'T', 'Tt', 'u_ax', 'u_th', 's']
            if variable not in variable_names:
                raise TypeError(
                    "The specified variable is not an output for the mixed out average. Use: 'P', 'Pt', 'T', 'Tt', 'u_ax', 'u_th' or 's'")
            index = variable_names.index(variable)
            for i in range(1, N_k + 1):
                mixedOut_result = self.get_mixed_out_average(j, [0, 1], [k_list[i - 1], k_list[i + 1]])
                avr.append(mixedOut_result[index])
        else:
            raise TypeError(
                "avr_type was not recognised. It must be one of the followings: 'areaAve', 'massFlowAve' or 'mixedOutAve'.")

        return np.array(avr)

    def get_area_average(self, variable: str, level: float, i_lim: list, k_lim: list) -> float:
        """ Method to compute the areaaverage of the variable "variable" in a plane defined by k=level·N_k and the
            limits set by 'i_lim' and 'j_lim'. It is based on the method  'pre_average' from this class.

            Parameters
            ----------
            variable : str
                Variable that is going to be averaged. It must be a key in the dictionary BlockCDF.flow_field.variables.
            level : float
                Number between 0 and 1 that determines the position of the grid plane
            i_lim : list
                List of 2 floats between 0 and 1 that determines the range i that is going to be used for the averaging.
            k_lim : list
                List of 2 floats between 0 and 1 that determines the range k that is going to be used for the averaging.

            Returns
            -------
            avr : float
                value of the area averaged quantity
            """

        [var, dA] = self.pre_average([variable], level, i_lim, k_lim)
        fun = var * dA
        avr = fun.sum() / dA.sum()
        return avr

    def get_mass_flow_average(self, variable: str, level: float, i_lim: list, k_lim: list) -> float:
        """ Method to compute the mass flow average of the variable "variable" in a plane defined by k=level·N_k and the
            limits set by 'i_lim' and 'j_lim'. It is based on the method  'pre_average' from this class.

            Parameters
            ----------
            variable : str
                Variable that is going to be averaged. It must be a key in the dictionary BlockCDF.flow_field.variables.
            level : float
                Number between 0 and 1 that determines the position of the grid plane
            i_lim : list
                List of 2 floats between 0 and 1 that determines the range i that is going to be used for the averaging.
            k_lim : list
                List of 2 floats between 0 and 1 that determines the range k that is going to be used for the averaging.

            Returns
            -------
            avr : float
                value of the area averaged quantity
            """

        [var, dA] = self.pre_average([variable, 'rho', 'Vx_stn'], level, i_lim, k_lim)
        rho = var[1]
        Vx = var[2]
        fun = var[0] * rho * Vx * dA
        d_m_dot = rho * Vx * dA
        avr = fun.sum() / d_m_dot.sum()
        return avr

    def get_mixed_out_average(self, level: float, i_lim: list, k_lim: list) -> list:
        """ Method to compute the mixed out average of the variable "variable" in a plane defined by k=level·N_k and the
            limits set by 'i_lim' and 'j_lim'. Everything is based on the stationary frame of reference. It is based on
            the method  'pre_average' from this class.

            The mixed out state is taken as that one that a viscous fluid will achieve after evolving through a constant
            section duct of inviscid walls. This means:
                - No viscous forces between streamlines (u_ax = uniform, u_r = 0, u_th = B·r; solid body rotation).
                - No heat transfer between streamlines (T = uniform).
                - Pressure gradient such that radial equilibrium is satisfied: dP/dr = rho·u_th**2/r.
                - Conservation of axial momentum: int(P + rho·u_ax**2)dA is conserved.
                - Conservation of angular momentum: int(rho·u_th·r·u_ax)dA is conserved.
                - Conservation of energy: int(Tt·rho·u_ax)dA is conserved.
                - Conservation of mass: int(rho·u_ax)dA is conserved.
            Note that the mixed out state is in general rotational and not uniform in energy nor entropy distributions.

            Parameters
            ----------
            level : float
                Number between 0 and 1 that determines the position of the grid plane
            i_lim : list
                List of 2 floats between 0 and 1 that determines the range i that is going to be used for the averaging.
            k_lim : list
                List of 2 floats between 0 and 1 that determines the range k that is going to be used for the averaging.

            Returns
            -------
            avr : list
                List containing the averaged values of the mixed out state. If they are not uniform it is the mass
                average.
                    - avr[0] = P [Pa]
                    - avr[1] = Pt [Pa]
                    - avr[2] = T [K]
                    - avr[3] = Tt [K]
                    - avr[4] = u_m [m/s]
                    - avr[5] = u_th [m/s] (evaluated at r_Euler = (r_tip+r_root)/sqrt(2) )
                    - avr[6] = s [J/K] (T_ref = 288.15 K, P_ref = 101325 Pa)
            """

        def function(y: np.ndarray) -> np.ndarray:
            """
            :parameter y: list
                List with the actual variables of the function:
                - y[0] = B
                - y[1] = u_ax
                - y[2] = T0
                - y[3] = P0
            :return res : list
                list of residuals of each of the equations:
                - res[0] = Angular momentum
                - res[1] = Axial momentum
                - res[2] = Mass flow
                - res[3] = Total energy
            """

            N = 100  # Number of radial divisions of the domain
            B = y[0]
            u_ax = y[1]
            T0 = y[2]
            P0 = y[3]
            r_fun = np.linspace(R_data[0], R_data[1], N)

            u_th = r_fun * B
            u_ax = u_ax * np.ones(N)
            T = T0 * np.ones(N)
            EXP = [np.exp(B ** 2 * (r_fun[i] ** 2 - r_fun[0] ** 2) / (2 * Rg * T0)) for i in range(N)]
            P = P0 * np.array(EXP)
            rho = P / (Rg * T)
            Tt = T + (u_th ** 2 + u_ax ** 2) / (2 * Cp)

            r_dual = np.zeros(N + 1)
            r_dual[0] = r_fun[0]
            r_dual[N] = r_fun[N - 1]
            r_dual[1:-1] = 0.5 * (r_fun[:-1] + r_fun[1:])
            dA_fun = (r_dual[1:] ** 2 - r_dual[:-1] ** 2) * 0.5 * R_data[2]

            # Angular momentum
            dm = rho * u_th * r_fun * u_ax * dA_fun
            m_fun = dm.sum()
            res_m = m_fun - m

            # Axial momentum:
            dp = (P + rho * u_ax ** 2) * dA_fun
            p_fun = dp.sum()
            res_p = p_fun / p - 1

            # Mass flow
            d_m_dot = rho * u_ax * dA_fun
            m_dot_fun = d_m_dot.sum()
            res_mdot = m_dot_fun / m_dot - 1

            # Total energy:
            dTt = Tt * d_m_dot
            Tt_fun = dTt.sum() / m_dot_fun
            res_Tt = Tt_fun / Tt_avr - 1

            return np.array([res_m, res_p, res_mdot, res_Tt])

        [var, dA] = self.pre_average(['rho', 'P', 'Tt_stn', 'Vx_stn', 'Vth_stn', 'r'], level, i_lim, k_lim)
        rho = var[0]
        P = var[1]
        Tt = var[2]
        Vx = var[3]
        Vth = var[4]
        r = var[5]

        # Axial momentum flux [N]:
        dp = (P + rho * Vx ** 2) * dA
        p = dp.sum()

        # Angular momentum flux [N·m]:
        dm = rho * Vth * r * Vx * dA
        m = dm.sum()

        # mass flux [Kg/s]:
        d_m_dot = rho * Vx * dA
        m_dot = d_m_dot.sum()

        # Tt mass average (total enthalpy flux) [J / ( J/(Kg·K) )]:
        dTt = Tt * d_m_dot
        Tt_avr = dTt.sum() / m_dot

        # Thermodynamic properties of the fluid:
        Cp = self.raw_data_grid[2][0]
        gamma = self.raw_data_grid[2][1]
        Cv = Cp / gamma
        Rg = Cp - Cv

        # Extract data about the constant (annular) section:
        R_data = [r[0, 0], r[0, r.shape[1] - 1],
                  dA.sum() / (0.5 * (r[0, r.shape[1] - 1] ** 2 - r[0, 0] ** 2))]  # r_hub, r_shroud, Delta_theta

        # Solve the non-linear system of equations:
        N_iter = 50  # max number of iterations
        tol = 0.0001  # Convergence tolerance
        factor = 0.5  # relaxation factor, helps stability
        y_norm = np.array([Vx[0, 0], Vx[0, 0], Tt_avr / 1.1, P[0, 0]])  # initial guess
        y0 = np.ones(4)
        delta = 0.0001  # Used to compute the derivative
        J = np.zeros([len(y0), len(y0)])
        for i in range(N_iter):
            f0 = function(y0 * y_norm)
            for j in range(len(y0)):
                y = deepcopy(y0)
                y[j] += delta
                J[:, j] = (function(y * y_norm) - f0) / delta

            y_new = y0 - factor * np.linalg.solve(J, f0)
            err = np.linalg.norm(y_new / y0 - 1)

            if err < tol:
                break
            elif i == N_iter - 1:
                raise TypeError("Newton-Raphson scheme is not converging")
            else:
                y0 = y_new

        # These are the parameters defining the mixed out state:
        y = y_new * y_norm
        B = y[0]
        u_m = y[1]
        T0 = y[2]
        P0 = y[3]

        # Compute averaged values:
        N = 50  # Number of radial positions to compute the mixed out state
        r_mixed_out = np.linspace(R_data[0], R_data[1], N)
        u_th = r_mixed_out * B
        u_ax = u_m * np.ones(N)
        T = T0 * np.ones(N)
        EXP = [B ** 2 * (r_mixed_out ** 2 - r_mixed_out[0] ** 2) / (2 * Rg * T0)]
        P = P0 * np.exp(EXP)
        rho = P / (Rg * T)
        Tt = T + (u_th ** 2 + u_ax ** 2) / (2 * Cp)
        Pt = P * (Tt / T) ** (gamma / (gamma - 1))

        r_dual = np.zeros(N + 1)
        r_dual[0] = r_mixed_out[0]
        r_dual[N] = r_mixed_out[N - 1]
        r_dual[1:-1] = 0.5 * (r_mixed_out[:-1] + r_mixed_out[1:])
        dA = (r_dual[1:] ** 2 - r_dual[:-1] ** 2) * 0.5 * R_data[2]

        dP = P * rho * u_ax * dA
        P = dP.sum() / m_dot

        dPt = Pt * rho * u_ax * dA
        Pt = dPt.sum() / m_dot

        r_Euler = np.sqrt((R_data[0] ** 2 + R_data[1] ** 2) / 2)
        u_th = B * r_Euler

        T_ref = 288.15
        P_ref = 101325
        s = Cp * np.log(T0 / T_ref) - Rg * np.log(P / P_ref)

        avr = [P, Pt, T0, Tt_avr, u_m, u_th, s]

        return avr

    def pre_average(self, variable: list, level: float, i_lim: list, k_lim: list) -> (list, np.ndarray):
        """ Method that supports all the other averaging methods. It takes the list of variable names 'variable: str',
            which must be keys in the "variable" attribute of BlockCFD.flow_field.
            This method extracts the grid differential areas and variables for a grid surfaces j = constant.

            Parameters
            ----------
            variable : list
                List of strings with the ariables to be extracted. They must be a key in the dictionary
                BlockCDF.flow_field.variables.
            level : float
                Number between 0 and 1 that determines the position of the grid plane
            i_lim : list
                List of 2 floats between 0 and 1 that determines the range i that is going to be used for the averaging.
            k_lim : list
                List of 2 floats between 0 and 1 that determines the range k that is going to be used for the averaging.

            Returns
            -------
            variables : list
                list of numpy nparrays containing the selected variables.
            dA : numpy.nparray
                grid of differential area associated with each grid point.
            """

        # Size the plane of averaging:
        N_i = self.N_i
        N_j = self.N_j
        N_k = self.N_k

        j = int(np.round((N_j - 1) * level))
        i_low = int(np.round((N_i - 1) * i_lim[0]))
        i_high = int(np.round((N_i - 1) * i_lim[1]))
        k_low = int(np.round((N_k - 1) * k_lim[0]))
        k_high = int(np.round((N_k - 1) * k_lim[1]))

        # Extract the variables at the section of interest:
        var = []
        for var_name in variable:
            try:
                var.append(self.flow_field.variables[var_name][i_low:i_high + 1, j, k_low:k_high + 1])
            except:
                raise TypeError("The variable specified does not exist, check this_object.flow_field.variables.keys")

        # Extract the coordinates of the plane
        rth = self.grid.rth[i_low:i_high + 1, j, k_low:k_high + 1]
        r = self.grid.r[i_low:i_high + 1, j, k_low:k_high + 1]

        # Compute the dA associated to every point of the plane
        [N_x, N_y] = rth.shape
        original_grid = [rth, r]
        dual_grid = [np.zeros([N_x + 1, N_y + 1]), np.zeros([N_x + 1, N_y + 1])]

        for k in range(2):
            # Corners:
            dual_grid[k][0, 0] = original_grid[k][0, 0]
            dual_grid[k][N_x, 0] = original_grid[k][N_x - 1, 0]
            dual_grid[k][0, N_y] = original_grid[k][0, N_y - 1]
            dual_grid[k][N_x, N_y] = original_grid[k][N_x - 1, N_y - 1]

            # Edges
            dual_grid[k][1:N_x, 0] = 0.5 * (original_grid[k][0:N_x - 1, 0] + original_grid[k][1:N_x, 0])
            dual_grid[k][1:N_x, N_y] = 0.5 * (original_grid[k][0:N_x - 1, N_y - 1] + original_grid[k][1:N_x, N_y - 1])
            dual_grid[k][0, 1:N_y] = 0.5 * (original_grid[k][0, 0:N_y - 1] + original_grid[k][0, 1:N_y])
            dual_grid[k][N_x, 1:N_y] = 0.5 * (original_grid[k][N_x - 1, 0:N_y - 1] + original_grid[k][N_x - 1, 1:N_y])

            # Middle of the domain:
            dual_grid[k][1:N_x, 1:N_y] = 0.25 * (
                        original_grid[k][0:N_x - 1, 0:N_y - 1] + original_grid[k][1:N_x, 0:N_y - 1] + original_grid[k][
                                                                                                      0:N_x - 1,
                                                                                                      1:N_y] +
                        original_grid[k][1:N_x, 1:N_y])

        # This was a test:
        # plt.figure()
        # plt.plot(original_grid[0],original_grid[1],'bo')
        # plt.plot(dual_grid[0],dual_grid[1],'r')
        # plt.plot(dual_grid[0].transpose(), dual_grid[1].transpose(), 'r')
        # plt.xlabel(r'$r\cdot\theta$')
        # plt.ylabel('r')
        # plt.title('Axis not to scale')

        dA = np.zeros([N_x, N_y])
        dA[:, :] = 0.5 * (
                    (dual_grid[0][1:, 1:] - dual_grid[0][:-1, 1:]) + (dual_grid[0][1:, :-1] - dual_grid[0][:-1, :-1])) * \
                   0.5 * ((dual_grid[1][1:, 1:] - dual_grid[1][1:, :-1]) + (
                    dual_grid[1][:-1, 1:] - dual_grid[1][:-1, :-1]))

        return var, dA
