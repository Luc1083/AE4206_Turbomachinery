import numpy as np
import scipy.optimize as opt
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


# Estimation of axial compressor performance along the blade height.
# Authors: ir. Evert Bunschoten, Dr. Andrea Giuffre', Dr. ir. Matteo Pini
# Reference: C. Freeman, N.A. Cumpsty - "A Method for the Prediction of Supersonic Compressor Blade Performance", 1989
# Delft University of Technology - All rights reserved


class CompressorEstimator:
    # Emulator class for quasi-1D compressor inlet for transonic conditions.

    _beta_blade:np.ndarray[float] = np.array([-60.4,-60.4,-60.4]) # Blade metal angles at root, mid, and tip [deg]
    _R_hub:float = 0.17 # Hub radius [m]
    _R_tip:float = 0.3  # Tip radius [m]
    _R_mean:float   # Mean radius [m]
    _Radius:np.ndarray[float]   # Radius array
    _beta_blade_interp:interp1d # Interpolator object for interpolation of the blade metal angles along blade span.

    _t_ratio:float = 0.05   # Thickness-to-chord ratio [-]
    _throat_area:np.ndarray[float]  # Array with throat area along the blade height.
    _t_ratio_blade:np.ndarray[float]    # Thickness-to-pitch ratio along the blade height.
    _t_ratio_interp:interp1d    # Interpolator object for interpolation of the ts ratio along blade height.

    _N_bl:int = 19  # Blade count [-]
    _val_pitch:np.ndarray[float]    # Blade pitch values along blade height.
    _soliditiy:float = 1.0  # Row solidity value applied for the whole blade.

    # Operating conditions
    _Pt:float = 276000  # Inlet stagnation pressure [Pa]
    _Tt:float = 340.0   # Inlet stagnation temperature [K]
    _mdot:float = 45.0  # Mass flow rate [kg s^-1]
    _rpm:float = 13177  # Rotor rotation speed [rpm]
    __started:bool=True # Started(True)/unstrated(False) condition.

    # Working fluid (air)
    _gamma:float=1.4    # Specific heat ratio
    _R_gas:float=287.05 # Gas constant [J kg^-1 K^-1]
    _cp:float = 1004.5  # Specific heat at constant pressure [J kg^-1 K^-1]

    # Local flow quantities
    _M_1:float      # Relative inflow Mach number.
    _beta_1:float   # Relative inflow angle [deg]
    _T_s1:float     # Static inflow temperature [K].
    _P_s1:float     # Static inflow pressure [Pa].

    # Quantities along blade height
    __Np_along_blade = 100  # Number of points to plot
    __incidence_along_blade = np.zeros(__Np_along_blade)    # Incidence angle along blade height [deg]
    __residual_along_blade = np.zeros(__Np_along_blade)     # Freeman control-volume equation residual along blade height.
    __M_1_along_blade = np.zeros(__Np_along_blade)          # Inflow Mach number distribution.
    __M_2_along_blade = np.zeros(__Np_along_blade)          # Throat Mach number distriubtion.

    

    def __init__(self):
        print("Defining a compressor emulator")
        self.__CompileGeom()

    def SetBetaBlade(self, val_beta:float, section:int=None):
        """Define the blade metal angle for a specific section of the 
        compressor blade or the whole blade.

        :param val_beta: blade metal angle in degrees.
        :type val_beta: float
        :param section: section at which to apply the given angle, defaults to None
        :type section: int, optional
        """
        if section == None:
            for i in range(3):
                self._beta_blade[i] = val_beta
        else:
            self._beta_blade[section] = val_beta
        self.__CompileGeom()

    def GetBetaBlade(self, section:int):
        """Return the blade metal angle value at a specified section.

        :param section: Blade section index (0:root, 1:mid, 2:tip)
        :type section: int
        :raises Exception: if section index is below 0 or above 2.
        :return: blade metal angle value at the specified section in degrees.
        :rtype: float
        """
        if section < 0 or section > 2:
            raise Exception("Section value should be 0, 1, or 2")
        return self._beta_blade[section]
    
    def SetHubRadius(self, val_hub_radius:float):
        self._R_hub = val_hub_radius
        self.__CompileGeom()

    def GetHubRadius(self):
        return self._R_hub
    
    def SetShroudRadius(self, val_shroud_radius:float):
        self._R_tip = val_shroud_radius
        self.__CompileGeom()

    def GetShroudRadius(self):
        return self._R_tip
    
    def SetThicknessRatio(self, val_t_ratio:float):
        self._t_ratio = val_t_ratio
        self.__CompileGeom()

    def GetThicknessRatio(self):
        return self._t_ratio
    
    def __CompileGeom(self):
        self._R_mean = 0.5*(self._R_tip + self._R_hub)
        self._Radius = np.array([self._R_hub, self._R_mean, self._R_tip])
        self._val_pitch = 2*np.pi*self._Radius / self._N_bl
        # Default chord length:10cm
        chord = 0.1
        t_LE = self._t_ratio * chord
        self._throat = (self._val_pitch * np.cos(np.deg2rad(self._beta_blade)))
        self._t_ratio_blade = t_LE / self._throat
        self._beta_blade_interp = interp1d(self._Radius, self._beta_blade,bounds_error=False)
        self._t_ratio_interp = interp1d(self._Radius, self._t_ratio_blade,bounds_error=False)
    
    def SetBladeCount(self, val_NB:int):
        self._N_bl = val_NB
        self.__CompileGeom()

    def GetBladeCount(self):
        return self._N_bl
    
    def ChangeOperatingCondition(self, val_Pt:float, val_Tt:float, val_mdot:float, val_rpm:float):
        self._Pt = val_Pt
        self._Tt = val_Tt 
        self._mdot = val_mdot 
        self._rpm = val_rpm 

    def SetPt(self, val_Pt:float):
        self._Pt = val_Pt
    def GetPt(self):
        return self._Pt 
    def SetTt(self, val_Tt:float):
        self._Tt = val_Tt
    def GetTt(self):
        return self._Tt
    def SetMdot(self, val_mdot:float):
        self._mdot = val_mdot
    def GetMdot(self):
        return self._mdot 
    def SetRPM(self, val_rpm:float):
        self._rpm = val_rpm 
    def GetRPM(self):
        return self._rpm
    
    def ChangeFluid(self, val_gamma:float, val_R:float):
        self._gamma = val_gamma 
        self._R_gas = val_R
        self._cp = self._gamma * self._R_gas / (self._gamma - 1)

    def GetBeta_1(self):
        return self._beta_1
    def GetMach_1(self):
        return self._M_1 
    
    def __compute_inlet_flow_conditions(self, R_loc):
        """
        Given the operating point specified in terms of mass flow rate and rotational speed and
        the inlet geometry specified in terms of blade height, mean radius and blade angle,
        compute the inlet relative Mach number and the flow angle.
        """
        mass_flow = self._mdot
        Area = 2*np.pi*R_loc * (self._R_tip - self._R_hub)
        omega = self._rpm * (2 * np.pi) / 60
        U = omega * R_loc
        Mach_abs = opt.fsolve(self.__massflow_isentropic, (0.5), args=(Area, mass_flow), full_output=False)[0]

        T = self._Tt / (1 + (self._gamma - 1) / 2 * Mach_abs ** 2)
        P = self._Pt / ((self._Tt / T) ** (self._gamma / (self._gamma - 1)))
        Vm = Mach_abs * np.sqrt(self._gamma * self._R_gas * T)
        Wm = Vm
        Wt = 0.0 - U
        W = np.sqrt(Wm ** 2 + Wt ** 2)
        Mach_rel = W / np.sqrt(self._gamma * self._R_gas * T)
        beta = np.arctan(Wt / Wm)

        self._M_1 = Mach_rel 
        self._beta_1 = np.rad2deg(beta)
        self._T_s1 = T 
        self._P_s1 = P 


    def __massflow_isentropic(self, p, *data):
        """
        Given the annulus area and the mass flow passing through it, compute the absolute Mach number, assuming
        isentropic flow and perfect gas.
        """
        A, m = data
        Mach = p

        m_computed = self._Pt * A / np.sqrt(self._R_gas * self._Tt) * Mach * np.sqrt(self._gamma) * (1 + (self._gamma - 1) / 2 * Mach ** 2) ** \
                    ((1 + self._gamma) / (2 * (1 - self._gamma)))
        res = (m_computed - m) / m

        return res
    
    def __freeman_cv(self, t_th, beta, beta_blade, Mach_rel_in, Mach_rel_out):
        """
        Compute left-hand and right-hand sides of the equation derived by Freeman and Cumpsty for the estimation
        of choking point in axial compressor cascades, assuming perfect gas.
        Inlet and outlet refer to the boundaries of the control volume, not of the cascade.
        """
        gamma = self._gamma
        lhs = ((1 + (gamma - 1) / 2 * (Mach_rel_out ** 2)) ** (- 1 / 2)) * \
            (1 + gamma * (Mach_rel_out ** 2) * (1 - t_th)) / (Mach_rel_out * (1 - t_th))
        rhs = ((1 + (gamma - 1) / 2 * (Mach_rel_in ** 2)) ** (- 1 / 2)) * \
            (np.cos(beta_blade) / np.cos(beta) + gamma * (Mach_rel_in ** 2) * np.cos(beta - beta_blade)) / Mach_rel_in
        
        return lhs, rhs   

    def __compute_post_shocks_flow_conditions(self, M_2, *data):
        """
        Compute Mach number at the outlet of the control volume defined by Freeman and Cumpsty.
        """
        t_th, beta, beta_blade, Mach_rel_in = data
        Mach_rel_out = M_2

        lhs, rhs = self.__freeman_cv(t_th, beta, beta_blade, Mach_rel_in, Mach_rel_out)
        res = (lhs - rhs) / lhs

        return res
    

    def ComputeLocalPerformance(self, radius_factor):
        val_radius = self._R_hub + radius_factor * (self._R_tip - self._R_hub)
        self.__compute_inlet_flow_conditions(val_radius)
        beta_blade_local = self._beta_blade_interp(val_radius)
        beta_blade_local_rad = np.deg2rad(beta_blade_local)
        t_th_local = self._t_ratio_interp(val_radius)

        
        M_1_local = self._M_1 
        beta_1_local = self._beta_1
        beta_1_local_rad = np.deg2rad(beta_1_local)

        val_indicence = np.abs(beta_1_local) - np.abs(beta_blade_local)
        #print(beta_1_local, beta_blade_local)
        theta = beta_1_local - beta_blade_local
        #print(theta)
        

        x0= 0.5*M_1_local
        result = opt.fsolve(self.__compute_post_shocks_flow_conditions, (x0),
                                args=(t_th_local, beta_1_local_rad, beta_blade_local_rad, M_1_local), full_output=True)

        residual = result[1]['fvec']
        M_2 = result[0][0]

        # if M_1_local > 1:
        #     if theta < 0:
        #         beta_shock = self.ComputeObiqueShock(M_1_local, -theta)[0]
        #     else:
        #         beta_shock = self.ComputeObiqueShock(M_1_local, theta)[0]
                
        #     beta_shock = np.rad2deg(beta_shock)
        #     if beta_shock < 90 and beta_shock > 0:
        #         print(beta_shock)

        return M_1_local, val_indicence, residual, M_2
    
    def ComputeObiqueShock(self, M_1, theta):
        beta_shock = opt.fsolve(self.__oblique_shock_relation, (np.deg2rad(theta)), args=(np.deg2rad(theta), M_1), full_output=False)
        return beta_shock 
    
    def __oblique_shock_relation(self, beta, theta, M_1):
        LHS = np.tan(theta)
        term_1 = 2 / np.tan(beta)
        term_2 = np.power(M_1*np.sin(beta), 2) - 1
        term_3 = np.power(M_1, 2) * (self._gamma + np.cos(2*beta))+2
        RHS = term_1*term_2/term_3
        return LHS - RHS
    
    def ComputePerformance(self):
        radius_factor_range = np.linspace(0, 1, self.__Np_along_blade)
        self.__incidence_along_blade = np.zeros(self.__Np_along_blade)
        self.__residual_along_blade = np.zeros(self.__Np_along_blade)
        self.__M_1_along_blade = np.zeros(self.__Np_along_blade)
        self.__M_2_along_blade = np.zeros(self.__Np_along_blade)
        
        for i_loc, val_radius in enumerate(radius_factor_range):
            M_1_local, val_incidence, residual, M_2 = self.ComputeLocalPerformance(val_radius)
            
            self.__incidence_along_blade[i_loc] = val_incidence
            self.__M_1_along_blade[i_loc] = M_1_local
            self.__residual_along_blade[i_loc] = residual
            self.__M_2_along_blade[i_loc] = M_2

    def PlotSpanWiseResults(self, ax_to_plot=None,current_location=None):
        if ax_to_plot == None:
            fig = plt.figure(figsize=[10,10])
            ax = plt.axes()
        else:
            ax = ax_to_plot

        spanwise_fraction_range = np.linspace(0, 1, self.__Np_along_blade)

        ax.plot(self.__M_1_along_blade, spanwise_fraction_range, 'm', label='Intake')
        ax.plot(self.__M_2_along_blade, spanwise_fraction_range, 'b', label='At throat')

        if current_location is not None:
            M1_interp = np.interp(current_location, np.linspace(0, 1, self.__Np_along_blade), self.__M_1_along_blade)
            M2_interp = np.interp(current_location, np.linspace(0, 1, self.__Np_along_blade), self.__M_2_along_blade)
            ax.plot(M1_interp, current_location, 'mo')
            ax.plot(M2_interp, current_location, 'bo')

        ax.set_ylabel(r"$\frac{r - r_h}{h}$",fontsize=20, rotation=0,labelpad=20)
        ax.set_yticks([])
        ax.set_xlabel(r"$M$",fontsize=20)
        ax.legend()
    def PlotBlade(self, blade_height_factor, ax_to_plot=None):
        val_radius = self._R_hub + blade_height_factor * ((self._R_tip - self._R_hub))
        val_beta_blade = self._beta_blade_interp(val_radius)

        x_range = np.linspace(0, 1 / self._soliditiy, 20)
        y_range_base = x_range * np.tan(np.deg2rad(-val_beta_blade))
        if ax_to_plot == None:
            fig = plt.figure()
            ax = plt.axes()
        else:
            ax = ax_to_plot 

        ax.plot(x_range, y_range_base, 'k')
        ax.plot(x_range, y_range_base+1, 'k')
        ax.plot(x_range, y_range_base-1, 'k')
        ax.set_xlabel(r"$\frac{c}{s}$", fontsize=20)
        ax.set_ylabel(r"$\frac{r\theta}{s}$", fontsize=20,rotation=0,labelpad=20)
        ax.set_xticks([])
        ax.set_yticks([])
        return val_beta_blade
    
    def PlotBladeAndFlow(self, blade_height_factor, ax_to_plot=None):
        M1, i, r_started, M2_started= self.ComputeLocalPerformance(blade_height_factor)

        beta_blade = self.PlotBlade(blade_height_factor=blade_height_factor, ax_to_plot=ax_to_plot)
        incidence = abs(self._beta_1) - abs(beta_blade)

        dX_Mach = M1*np.cos(np.deg2rad(self._beta_1))
        dY_Mach = M1*np.sin(np.deg2rad(self._beta_1))

        dX_Mach2 = M2_started*np.cos(np.deg2rad(-beta_blade))
        dY_Mach2 = M2_started*np.sin(np.deg2rad(-beta_blade))

        if ax_to_plot == None:
            fig = plt.figure()
            ax = plt.axes()
        else:
            ax = ax_to_plot
        title_str = ""
        if abs(r_started) > 1e-6:
            title_str += "Chocked!"
            #ax.set_title("Chocked!")
        elif incidence > 10:
            title_str += "High incidence!"
        
        ax.set_title(title_str)
        ax.arrow(x=0, y=0.5, dx=dX_Mach2, dy=dY_Mach2,color='b',width=0.04)
        ax.arrow(x=-dX_Mach, y=0.5+dY_Mach, dx=dX_Mach, dy=-dY_Mach,color='m',width=0.04)
        ax.text(x=-dX_Mach, y=0.5+dY_Mach, s=(r"$M_1=%.2f$"%M1),horizontalalignment='center',verticalalignment='top',color='m')
        ax.text(x=dX_Mach2, y=0.6+dY_Mach2, s=(r"$M_2=%.2f$"%M2_started),horizontalalignment='center',verticalalignment='bottom',color='b')
        
        ax.set_xlim([-2*dX_Mach, 1.2])
        
        ax.set_aspect('equal')


