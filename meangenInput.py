import time
import warnings
from sys import platform
import os
import subprocess
import meanline as ml
import numpy as np

import matplotlib.pyplot as plt


def angle_between_vectors(v1, v2=np.array([1, 0])):
    return np.arccos(v1.dot(v2) / np.linalg.norm(v1) / np.linalg.norm(v2))


class MeangenCompressorInput:
    def __init__(self, fan: ml.Fan, force_axial_chords: bool = False, force_axial_gaps: bool = False,
                 force_blockage_factors: bool = False, force_deviation: bool = True, force_eta: bool = False,
                 force_incidence: bool = True, force_twist_factor: bool = True, force_blade_tc_loc: bool = True,
                 force_q0: bool = True, exec_extension=""):
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
        :param exec_extension: file extension for executable
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
        self.force_Q0 = force_q0
        self.exec_extension = exec_extension

    def generate_input_file(self):
        with open(f"{os.getcwd()}/meangen.in", "w") as temp:
            temp.write("C\n")  # compressor
            temp.write("AXI\n")  # axial
            temp.write(f"{self.fan.R_air :.3f}     {self.fan.gamma :.3f}\n")  # sets R and Gamma
            temp.write(f"{self.fan.P0_inlet / 1e5 :.3f} {self.fan.T0_inlet :.3f}\n")  # provide inlet total values
            temp.write("1\n")  # set number of stages
            temp.write("M\n")  # specify design at mean line
            temp.write(f"{self.fan.omega}\n")  # fan RPM
            temp.write(f"{self.fan.mdot :.3f}\n")  # fan mass flow
            temp.write("A\n")  # specify that duty coefficients will be provided

            temp.write(
                f"{self.fan.R_mean :.3f}, {self.fan.theta_rotor_distribution[self.fan.rotor_mean_idx]:.3f}, "
                f"{self.fan.psi_mean :.3f}\n")  # provide duty coeffs: reaction, flow, loading
            temp.write("A\n")  # set method to specify tip radius
            temp.write(f"{self.fan.r_tip :.5}\n")  # specify tip radius in meters

            # set our own chord lengths or use meangen values
            if self.force_axial_chords:
                temp.write(f"{self.fan.c_mean_rotor :.5} {self.fan.c_mean_stator :.5}\n")
            else:
                temp.write("       0.050       0.040\n")

            # set our own row/stage gaps or use defaults
            # NOTE: I don't know if we compute this as of right now
            if self.force_axial_gaps:
                raise NotImplementedError("couldn't find it lol")
            else:
                temp.write("       0.250       0.500\n")

            # set our own blockage factors or use defaults. Factors provided at LE of row 1, and TE of row 2.
            # NOTE: couldn't find it
            if self.force_blockage_factors:
                raise NotImplementedError("couldn't find it lol")
            else:
                temp.write("   0.00000   0.00000     \n")

            # set our own stage efficiency or use default
            if self.force_eta:
                temp.write(f"{self.fan.eta_tt_estimated}\n")
                # temp.write(f"    {self.fan.eta}")
            else:
                temp.write("       0.900\n")

            # set our own deviation angles or use default
            if self.force_deviation:
                temp.write(f"{self.fan.delta_rotor[self.fan.rotor_mean_idx] :.3f} "  # do not remove space
                                f"{self.fan.delta_stator[self.fan.stator_mean_idx] :.3f}\n")
            else:
                temp.write("   5.000   5.000\n")

            # set our own incidence angles or use default
            if self.force_incidence:
                temp.write(f"{self.fan.i_rotor[self.fan.rotor_mean_idx] :.3f} "
                                f"{self.fan.i_stator[self.fan.stator_mean_idx] :.3f}\n")
            else:
                temp.write("  -2.000  -2.000\n")

            # twist fraction, 1 is free vortex, 0 is prismatic
            if self.force_twist_fact:
                temp.write(f"{self.fan.n :.4f}\n")
            else:
                temp.write("   1.00000      \n")

            # rotation option
            temp.write("N\n")  # We're not rotating sections rn I think, so I'm tactically ignoring this.

            if self.force_Q0:
                # Rotor angles

                rot_angle_in = np.rad2deg(
                    angle_between_vectors(np.array([(self.fan.c_hub_rotor - self.fan.c_tip_rotor) / 2,
                                                    self.fan.r_tip - self.fan.r_hub_inlet_rotor])))
                rot_angle_out = 180 - rot_angle_in
                # stator angles
                stat_angle_in = np.rad2deg(
                    angle_between_vectors(np.array([(self.fan.c_hub_stator - self.fan.c_tip_stator) / 2,
                                                    self.fan.r_tip - self.fan.r_hub_inlet_rotor])))
                stat_angle_out = 180 - stat_angle_in

                temp.write(f" {rot_angle_in :.3f} {rot_angle_out :.3f}\n")
                temp.write(f" {stat_angle_in :.3f} {stat_angle_out :.3f}\n")

            else:
                temp.write(
                    "  88.000  92.000\n")  # accepts default "Q0" angles for blade row 1. I have no idea what that is.
                temp.write(
                    "  92.000  88.000\n")  # accepts default "Q0" angles for blade row 2. I have no idea what that is.
            temp.write("N\n")  # again tells the thing to not rotate sections.
            temp.write("Y\n")  # output to stagen.dat

            if self.force_blade_tc_loc:
                warnings.warn("I don't think we specify the location of the max T/C, "
                              "replace the 0.4 below if I'm wrong\n")
                temp.write("N\n")  # do not accept defaults
                # rotor
                temp.write(f"{self.fan.t_c_rotor[0]} {0.4 :.5}\n")  # root
                temp.write(f"{self.fan.t_c_rotor[self.fan.rotor_mean_idx]} {0.4 :.5}\n")  # mid
                temp.write(f"{self.fan.t_c_rotor[-1]} {0.4 :.5}\n")  # tip

                temp.write("N\n")  # do not accept defaults
                temp.write(f"{self.fan.t_c_stator[0] :.5} {0.45 :.5}\n")  # root
                temp.write(f"{self.fan.t_c_stator[self.fan.stator_mean_idx] :.5} {0.45 :.5}\n")  # mid
                temp.write(f"{self.fan.t_c_stator[-1] :.5} {0.45 :.5}\n")  # tip
                warnings.warn("I'm assuming the sections 1, 2, and 3 map to the root, mean, and tip sections. "
                              "If you see this warning that means I have not checked it yet.\n")
            else:
                temp.write("N\n")  # do not accept defaults
                # rotor
                temp.write(f"  0.0750  0.4000\n")  # root
                temp.write(f"  0.0750  0.4000\n")  # mid
                temp.write(f"  0.0750  0.4000\n")  # tip

                temp.write("N\n")  # do not accept defaults
                temp.write("  0.1000  0.4500\n")  # root
                temp.write("  0.1000  0.4500\n")  # mid
                temp.write("  0.1000  0.4500\n")  # tip

        print("Meangen input file written")

    def run_meangen(self):
        p = subprocess.Popen(f"{os.getcwd()}/execs/meangen-17.4{self.exec_extension}", stdin=subprocess.PIPE, shell=True)
        p.stdin.write(bytes("F", "utf-8"))
        p.stdin.flush()
        p.stdin.close()
        p.wait()
        if p.returncode != 0:
            raise ChildProcessError("Goofy ahh Meangen code died")
        else:
            print("*********\nMeangen finished\n*********")


class RunCFD:
    def __init__(self, fan: ml.Fan, result_dir: str, overwrite=True):
        if platform == "linux":
            self.exec_extension = ""
        elif platform == "win32" or platform == "win64":
            self.exec_extension = ".exe"
        else:
            raise OSError("Yeah I dunno what kinda OS you're running but it's wrong :)")

        self.meangen_inp = MeangenCompressorInput(fan, exec_extension=self.exec_extension)
        self.result_dir = os.getcwd() + "/" + result_dir

        if os.path.isdir(self.result_dir):
            print("Overwriting existing result dir, it contains:")
            print(os.listdir(self.result_dir))
            if not overwrite:
                warnings.warn("Code is not set to overwrite so we're done")
                exit()
            os.rmdir(self.result_dir)
        os.mkdir(self.result_dir)

    def run_all(self):
        self.generate_meangen_input()
        self.run_stagen()
        self.run_multall()
        self.post_process()

    def generate_meangen_input(self):
        force_vals = [i for i in self.meangen_inp.__dict__ if "force" in i]
        for fv in force_vals:
            print(f"{fv} is set to {self.meangen_inp.__dict__[fv]}")
        print("Make sure you're happy with this")
        time.sleep(2)
        self.meangen_inp.generate_input_file()
        self.meangen_inp.run_meangen()

    def run_stagen(self):
        # figure out what the fuck
        p = subprocess.Popen(f"{os.getcwd()}/execs/stagen-18.1{self.exec_extension}", stdin=subprocess.PIPE, shell=True)
        p.stdin.write(bytes("Y", "utf-8"))
        p.stdin.flush()
        p.stdin.close()
        p.wait()
        if p.returncode != 0:
            raise ChildProcessError("Goofy ahh Stagen code died")
        else:
            print("*********\nStagen finished\n*********")

    def run_multall(self):
        if not os.path.isfile("intype"):
            with open(f"{os.getcwd()}/intype", "w") as f:
                f.write("N")
        p = subprocess.Popen(f"{os.getcwd()}/execs/multall-open-20.9{self.exec_extension} <stage_new.dat >results.out", shell=True)

    def post_process(self):
        # move some junk as well
        raise NotImplementedError()

    def refine_mesh(self):
        # change IM and KM in stagen.dat
        pass


if __name__ == "__main__":
    f = ml.Fan(Mach_inlet=0.6, AR_rotor=5, AR_stator=3, taper_rotor=3, taper_stator=0.802, n=0.72, no_blades_rotor=40,
                 no_blades_stator=38, beta_tt=1.6, P0_cruise=39513.14, T0_cruise=250.13, mdot=80, omega=5000,
                 hub_tip_ratio=0.5, gamma=1.4, R_air=287, eta_tt_estimated=0.9, row_chord_spacing_ratio=0.5,
                 lieblein_model=ml.Lieblein_Model(),
                 profile="NACA-65", methodology="free vortex")
    cfd = RunCFD(f, "run")
    # cfd.meangen_inp.force_Q0 = 0
    cfd.run_all()
