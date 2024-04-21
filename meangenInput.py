import warnings
from sys import platform
import os
import subprocess
import meanline as ml


class MeangenCompressorInput:
    def __init__(self, fan: ml.Fan, force_axial_chords: bool = False, force_axial_gaps: bool = False,
                 force_blockage_factors: bool = False, force_deviation: bool = True, force_eta: bool = False,
                 force_incidence: bool = True, force_twist_factor: bool = True, force_blade_tc_loc: bool = True):
        """
        Constructs Meangen input file. Requires fan object as input. Other arguments are optional, and determine whether
        the values from the fan object are forced or meangen defaults are used.
        :param fan: Fan object to analyse
        :param force_axial_chords: Bool, default False as it is not implemented
        :param force_axial_gaps: Bool, default False as it is not implemented
        :param force_blockage_factors: Bool, default True
        :param force_deviation: Bool, default True
        :param force_eta: Bool, default False as it is not implemented
        :param force_incidence: Bool, default True
        :param force_twist_factor: Bool, default True
        :param force_blade_tc_loc: Bool, default True. Only partially implemented
        """
        self.force_blade_tc_loc = force_blade_tc_loc
        self.force_twist_fact = force_twist_factor
        self.force_incidence = force_incidence
        self.force_blockage_factors = force_blockage_factors
        self.force_deviation = force_deviation
        self.force_eta = force_eta
        self.force_axial_gaps = force_axial_gaps
        self.fan = fan
        self.force_axial_chords = force_axial_chords

        if platform == "linux":
            self.str_add = ""
        elif platform == "win32" or platform == "win64":
            self.str_add = ".exe"
        else:
            raise OSError("Yeah I dunno what kinda OS you're running but it's wrong :)")

    def generate_input_file(self):
        with open("./meangen.in", "w") as self.temp:
            self.temp.write("C\n")  # compressor
            self.temp.write("AXI\n")  # axial
            self.temp.write(f"{self.fan.R_air :.3f}     {self.fan.gamma :.3f}\n")  # sets R and Gamma
            self.temp.write(f"{self.fan.P0_inlet / 1e5 :.3f} {self.fan.T0_inlet :.3f}\n")  # provide inlet total values
            self.temp.write("1\n")  # set number of stages
            self.temp.write("M\n")  # specify design at mean line
            self.temp.write(f"{self.fan.omega}\n")  # fan RPM
            self.temp.write(f"{self.fan.mdot :.3f}\n")  # fan mass flow
            self.temp.write("A\n")  # specify that duty coefficients will be provided

            self.temp.write(
                f"{self.fan.R_mean :.3f}, {self.fan.theta_rotor_distribution[self.fan.rotor_mean_idx]:.3f}, "
                f"{self.fan.psi_mean :.3f}\n")  # provide duty coeffs: reaction, flow, loading
            self.temp.write("A\n")  # set method to specify tip radius
            self.temp.write(f"{self.fan.r_tip :.5}\n")  # specify tip radius in meters

            # set our own chord lengths or use meangen values
            if self.force_axial_chords:
                self.temp.write(f"{self.fan.c_mean_rotor :.5} {self.fan.c_mean_stator :.5}\n")
            else:
                self.temp.write("       0.050       0.040\n")

            # set our own row/stage gaps or use defaults
            # NOTE: I don't know if we compute this as of right now
            if self.force_axial_gaps:
                raise NotImplementedError("couldn't find it lol")
            else:
                self.temp.write("       0.250       0.500\n")

            # set our own blockage factors or use defaults. Factors provided at LE of row 1, and TE of row 2.
            # NOTE: couldn't find it
            if self.force_blockage_factors:
                raise NotImplementedError("couldn't find it lol")
            else:
                self.temp.write("   0.00000   0.00000     \n")

            # set our own stage efficiency or use default
            if self.force_eta:
                raise NotImplementedError("Could not find it lol")
                # self.temp.write(f"    {self.fan.eta}")
            else:
                self.temp.write("       0.900\n")

            # set our own deviation angles or use default
            if self.force_deviation:
                self.temp.write(f"{self.fan.delta_rotor[self.fan.rotor_mean_idx] :.3f} "  # do not remove space
                                f"{self.fan.delta_stator[self.fan.stator_mean_idx] :.3f}\n")
            else:
                self.temp.write("   5.000   5.000\n")

            # set our own incidence angles or use default
            if self.force_incidence:
                self.temp.write(f"{self.fan.i_rotor[self.fan.rotor_mean_idx] :.3f} "
                                f"{self.fan.i_stator[self.fan.stator_mean_idx] :.3f}\n")
            else:
                self.temp.write("  -2.000  -2.000\n")

            # twist fraction, 1 is free vortex, 0 is prismatic
            if self.force_twist_fact:
                self.temp.write(f"{self.fan.n :.4f}\n")
            else:
                self.temp.write("   1.00000      \n")

            # rotation option (we don't)
            self.temp.write("N\n")  # We're not rotating sections rn I think, so I'm tactically ignoring this.

            self.temp.write(
                "  88.000  92.000\n")  # accepts default "Q0" angles for blade row 1. I have no idea what that is.
            self.temp.write(
                "  92.000  88.000\n")  # accepts default "Q0" angles for blade row 2. I have no idea what that is.
            self.temp.write("N\n")  # again tells the thing to not rotate sections.
            self.temp.write("Y\n")  # output to stagen.dat

            if self.force_blade_tc_loc:
                warnings.warn("I don't think we specify the location of the max T/C, "
                              "replace the 0.4 below if I'm wrong\n")
                self.temp.write("N\n")  # do not accept defaults
                # rotor
                self.temp.write(f"{self.fan.t_c_rotor[0]} {0.4 :.5}\n")  # root
                self.temp.write(f"{self.fan.t_c_rotor[self.fan.rotor_mean_idx]} {0.4 :.5}\n")  # mid
                self.temp.write(f"{self.fan.t_c_rotor[-1]} {0.4 :.5}\n")  # tip

                self.temp.write("N\n")  # do not accept defaults
                self.temp.write(f"{self.fan.t_c_stator[0] :.5} {0.45 :.5}\n")  # root
                self.temp.write(f"{self.fan.t_c_stator[self.fan.stator_mean_idx] :.5} {0.45 :.5}\n")  # mid
                self.temp.write(f"{self.fan.t_c_stator[-1] :.5} {0.45 :.5}\n")  # tip
                warnings.warn("I'm assuming the sections 1, 2, and 3 map to the root, mean, and tip sections. "
                              "If you see this warning that means I have not checked it yet.\n")
            else:
                self.temp.write("N\n")  # do not accept defaults
                # rotor
                self.temp.write(f"  0.0750  0.4000\n")  # root
                self.temp.write(f"  0.0750  0.4000\n")  # mid
                self.temp.write(f"  0.0750  0.4000\n")  # tip

                self.temp.write("N\n")  # do not accept defaults
                self.temp.write("  0.1000  0.4500\n")  # root
                self.temp.write("  0.1000  0.4500\n")  # mid
                self.temp.write("  0.1000  0.4500\n")  # tip

        print("Meangen input file written")

    def run_meangen(self):
        cmd = f"./execs/meangen-17.4{self.str_add}"
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, shell=True)
        p.stdin.write(bytes("F", "utf-8"))
        p.stdin.flush()
        p.stdin.close()
        p.wait()
        if p.returncode != 0:
            raise ChildProcessError("Goofy ahh Meangen code died")
        else:
            print("Meangen finished, moving files")


class RunCFD:
    def __init__(self):
        ...

    def run_stagen(self):
        ...


if __name__ == "__main__":
    f = ml.Fan(Mach_inlet=0.6, AR_rotor=10, AR_stator=10, taper_rotor=0.7, taper_stator=0.5, n=1, no_blades_rotor=30,
               no_blades_stator=60, beta_tt=1.6, P0_cruise=39513.14, T0_cruise=250.13, mdot=80, omega=5000,
               hub_tip_ratio=0.3, gamma=1.4, R_air=287, eta_tt_estimated=0.9, Cp_air=1006, Cv_air=715.9,
               row_chord_spacing_ratio=0.5, lieblein_model=ml.Lieblein_Model(), profile="NACA-65",
               methodology="controlled vortex", rho=1, dyn_visc=1e-6)
    mci = MeangenCompressorInput(f)
    mci.generate_input_file()
    mci.run_meangen()
