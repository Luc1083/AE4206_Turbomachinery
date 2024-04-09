import numpy as np
import matplotlib.pyplot as plt
from Components._002_blade_row import BladeRow
from Components.full_machine_plotter import PlotMachine


class Turbomachine:
    def __init__(self, grid_file='grid_out', flow_file='flow_out'):
        self.grid_name = grid_file  # TODO: Make this impossible to change after instantiating, the user should create another object
        self.flow_name = flow_file

        self.raw_data_grid = self.get_grid_data()
        N_i = self.raw_data_grid[1][0]
        N_j = self.raw_data_grid[1][1]
        N_k = self.raw_data_grid[1][2]
        self.raw_data_flow = self.get_flow_data(N_i,N_j,N_k)

        [self.N_rows, j_mix, ref_values] = self.characterize_machine()

        [data_sets_grid, data_sets_flow] = self.split_stages(j_mix)

        self.rows = []
        for i in range(self.N_rows):
            self.rows.append(BladeRow(data_sets_grid[i], data_sets_flow[i], ref_values))

        self.N_iter = self.raw_data_flow[0][0] # Number of the iteration where the flow field is known
        self.plot = PlotMachine(self)

    def get_grid_data(self):
        print('Reading Fortran binary raw data from "' + self.grid_name + '"...')
        f = open(self.grid_name, 'rb')

        data = []
        # This is the expected type of raw_data_grid: 0='int32' 1='float32' (Done by looking at Mulltall-open-20.9.f source code)
        cheat_list = [0, 0, 1, 0, 1, 0, 1]
        #   0. Max iterations number
        #   1. Dimensions of the domain
        #   2. Fluid properties (Cp & gamma)
        #   3. 'INDLETE' (Ind_LE_TE: 1 if it is LE or TE, 0 if not)
        #   4. W_rad (Rotational speed in rad/s)
        #   5. N_blade (Number of blades in the row)
        #   6. coordinates (several arrays, actually i>7)

        i = 0
        while True:
            n_val = np.fromfile(f, dtype='int32', count=1).tolist()
            if len(n_val) == 0:
                break
            else:
                n_val = n_val[0]

            if cheat_list[i] == 0:
                val = np.fromfile(f, dtype='int32', count=round(n_val/4)).tolist()
            else:
                val = np.fromfile(f, dtype='float32', count=round(n_val/4)).tolist()
            n_val_check = np.fromfile(f, dtype='int32', count=1).tolist()[0]

            if i < 6:
                i += 1

            if not n_val == n_val_check:
                print('n_val: ' + str(n_val))
                print('n_val_check: ' + str(n_val_check))
                raise Exception("The file is not behaving as expected. The code is not valid.")

            data.append(val)

        f.close()
        print('Done, "' + self.grid_name + '" closed.')
        return data

    def get_flow_data(self,N_i,N_j,N_k):
        print('Reading Fortran binary raw data from "' + self.flow_name + '"...')
        f = open(self.flow_name, 'rb')

        data = []
        # Expected structure (Based on Multall-open-20.9.f source code)
        #   0. N_step (iteration when printing results)
        #   1. rho(x,y,z): written as a line vector iterating i the quickest and k the slowest (all variables like this)
        #   2. rho·Vx (x,y,z)
        #   3. rho·Vr (x,y,z)
        #   4. rho·Vth (x,y,z)
        #   5. rho·e (x,y,z): Internal energy times density
        #   6. ROSUB in Multall-open-20.9.f It looks like a correction for non-ideal fluids. It is ignored in convert-to-tecplot.f
        #   7. Q in convert-to-tecplot.f and it is ignored. The information here depends on the flow model
        #   8. QQ in convert-to-tecplot.f and it is ignored. The information here depends on the flow model

        while True:
            n_val = np.fromfile(f, dtype='int32', count=1).tolist()
            if len(n_val) == 0:
                break
            else:
                n_val = n_val[0]

            if round(n_val/4) == 1:
                val = np.fromfile(f, dtype='int32', count=round(n_val / 4)).tolist()
            else:
                val = np.fromfile(f, dtype='float32', count=round(n_val / 4)).tolist()
            n_val_check = np.fromfile(f, dtype='int32', count=1).tolist()[0]

            if not n_val == n_val_check:
                print('n_val: ' + str(n_val))
                print('n_val_check: ' + str(n_val_check))
                raise Exception("The file is not behaving as expected. The code is not valid.")

            data.append(val)

        f.close()
        print('Done, "' + self.flow_name + '" closed.')
        return data

    def characterize_machine(self):
        N_j = self.raw_data_grid[1][1]
        omega = self.raw_data_grid[4]

        omega_test = omega[0]
        J_mix = [0]  # J coordinate of the interface between blade rows (j is first one, j+1 is the next one)
        for j in range(N_j):
            if not (omega[j] == omega_test):
                omega_test = omega[j]
                J_mix.append(j)

        J_mix.append(N_j)  # This is pointing at a point out of the domain (Phython counts from 0 and Fortran from 1,
        # N_j is regarding Fortran) but this way I can take my stages from j=j_mix(k) until
        # j=j_mix(k+1) (both included) and don't miss any points.

        N_rows = len(J_mix) - 1  # Number of rows in the machine

        # Obtain inlet stagnation pressure and temperature (reference values)
        Cp = self.raw_data_grid[2][0]
        gamma = self.raw_data_grid[2][1]
        Cv = Cp / gamma
        Rg = Cp - Cv

        rho = self.raw_data_flow[1][0]
        rhoVx = self.raw_data_flow[2][0]
        rhoVr = self.raw_data_flow[3][0]
        rhoVth = self.raw_data_flow[4][0]
        rhoE = self.raw_data_flow[5][0]

        Vx = rhoVx / rho
        Vr = rhoVr / rho
        Vth = rhoVth / rho
        e = rhoE / rho

        V = np.sqrt(Vx ** 2 + Vr ** 2 + Vth ** 2)
        T = (e - V ** 2 / 2) / Cv
        P = Rg * T * rho
        Tt = T + V ** 2 / (2 * Cp)
        Pt = P * (Tt / T) ** (gamma / (gamma - 1))

        reference_values = [Tt, Pt]

        return N_rows, J_mix, reference_values

    def split_stages(self, j_mix):
        # Split the grid raw data in different raw data sets cutting between j_mix and j_mix+1
        # Grid
        data_sets_grid = []
        for j in range(self.N_rows):
            cutted_data = [self.raw_data_grid[0]]
            cutted_data.append([self.raw_data_grid[1][0], j_mix[j + 1] - j_mix[j], self.raw_data_grid[1][2]])
            cutted_data.append(self.raw_data_grid[2])
            cutted_data.append(self.raw_data_grid[3][j_mix[j]: j_mix[j + 1]])
            cutted_data.append(self.raw_data_grid[4][j_mix[j]: j_mix[j + 1]])
            cutted_data.append(self.raw_data_grid[5][j_mix[j]: j_mix[j + 1]])

            N_k = self.raw_data_grid[1][2]
            aux_coord = self.raw_data_grid[6 + j_mix[j] * N_k: 6 + j_mix[j + 1] * N_k]
            for i in range(len(aux_coord)):
                cutted_data.append(aux_coord[i])

            data_sets_grid.append(cutted_data)

        # Flow field
        N_i = self.raw_data_grid[1][0]
        N_j = self.raw_data_grid[1][1]
        N_k = self.raw_data_grid[1][2]
        data_sets_flow = []

        for j in range(self.N_rows):
            cutted_data = [self.raw_data_flow[0]]
            for i in range(1,len(self.raw_data_flow)):
                aux_property = []
                for k in range(N_k):
                    # The data is stored in a vector such that the i,k,k correspondence is index=i+j·Ni+k·(Ni·Nj), with the indices starting at 0 and finishing at N-1
                    aux_property.append( self.raw_data_flow[i][0 + j_mix[j]*N_i + k*(N_i*N_j) : 0 + j_mix[j+1]*N_i + k*(N_i*N_j)])
                aux_property = [item for sublist in aux_property for item in sublist]
                cutted_data.append( aux_property )

            data_sets_flow.append(cutted_data)

        return data_sets_grid, data_sets_flow

    def gen_ParaView_input(self):
        self.gen_passages_ParaView()
        self.gen_blades_ParaView()

    def gen_blades_ParaView(self):
        print('Generating extra blade geometry and printing to .dat file, this might take couple minutes...')
        N_rows = self.N_rows
        print('  ' + str(N_rows) + ' blade row(s) found.')

        # Ensure that the geometry is loaded:
        for i in range(N_rows):
            if not len(self.rows[i].passage_instances) == self.rows[i].N_instances:
                print('Updating geometry in row ' + str(i) + ', this might take some time...')
                self.rows[i].generate_extra_geometry()

        file = open('ParaView_TecPlotInterpreter_blades.dat', 'wt')

        # Write ParaView file title
        file.write(' TITLE = "Run generated by PostPy with ' + self.grid_name + ' and ' + self.flow_name + ' (only blade data)"\n')

        # Get the names of the variable3s that are going to be plotted
        variables = self.rows[0].blade_instances[0].flow_field.variables
        computed_variables_names = list(variables.keys())
        N_var = len(computed_variables_names)
        variable_names = 'Variables = "X","Y","Z"'
        for i in range(N_var):
            variable_names += ',"' + computed_variables_names[i] + '"'
        variable_names += '\n'
        file.write(variable_names)

        # Write each block
        for i in range(N_rows):
            N_instances = len(self.rows[i].blade_instances)
            print('  Blade row ' + str(i) + ': printing ' + str(N_instances) + ' blades.')

            N_i = self.rows[i].blade_original.N_i
            N_j = self.rows[i].blade_original.N_j
            N_k = self.rows[i].blade_original.N_k

            for j in range(N_instances):
                print('   writing blade ' + str(j+1) + '/' + str(N_instances) + '...')
                # Write header of the block:
                header = 'ZONE T = "Row ' + str(i) + ', blade ' + str(j) + '" I= ' + str(N_i) + ' J= ' + str(N_j) + ' K= ' + str(N_k) + '\n'

                # Create the block of data we want to write:
                data = np.zeros([N_i * N_j * N_k, N_var + 3])
                variables = self.rows[i].blade_instances[j].flow_field.variables

                data[:, 0] = np.reshape(self.rows[i].blade_instances[j].grid.x, -1, order='F')  # x
                data[:, 1] = np.reshape(self.rows[i].blade_instances[j].grid.y, -1, order='F')  # y
                data[:, 2] = np.reshape(self.rows[i].blade_instances[j].grid.z, -1, order='F')  # z

                for k in range(N_var):
                    data[:,3+k] = np.reshape(variables[computed_variables_names[k]], -1, order='F')

                # Write down the chink of data:
                np.savetxt(fname=file, X=data, fmt='%10.5e', header=header, comments='')

        file.close()
        print('Done! your file is "ParaView_TecPlotInterpreter_blades.dat"')
        print(' ')

    def gen_passages_ParaView(self):
        print('Generating extra passage geometry and printing to .dat file, this might take couple minutes...')
        N_rows = self.N_rows
        print('  ' + str(N_rows) + ' blade row(s) found.')

        # Ensure that the geometry is loaded:
        for i in range(N_rows):
            if not len(self.rows[i].passage_instances) == self.rows[i].N_instances:
                print('Updating geometry in row ' + str(i) + ', this might take some time...')
                self.rows[i].generate_extra_geometry()

        file = open('ParaView_TecPlotInterpreter_passages.dat', 'wt')

        # Write ParaView file title
        file.write(
            ' TITLE = "Run generated by PostPy with ' + self.grid_name + ' and ' + self.flow_name + ' (only blade data)"\n')

        # Get the names of the variable3s that are going to be plotted
        variables = self.rows[0].passage_instances[0].flow_field.variables
        computed_variables_names = list(variables.keys())
        N_var = len(computed_variables_names)
        variable_names = 'Variables = "X","Y","Z"'
        for i in range(N_var):
            variable_names += ',"' + computed_variables_names[i] + '"'
        variable_names += '\n'
        file.write(variable_names)

        # Write each block
        for i in range(N_rows):
            N_instances = len(self.rows[i].passage_instances)
            print('  Blade row ' + str(i) + ': printing ' + str(N_instances) + ' passage(s).')

            N_i = self.rows[i].passage_original.N_i
            N_j = self.rows[i].passage_original.N_j
            N_k = self.rows[i].passage_original.N_k

            for j in range(N_instances):
                print('   writing passage ' + str(j + 1) + '/' + str(N_instances) + '...')
                # Write header of the block:
                header = 'ZONE T = "Row ' + str(i) + ', passage ' + str(j) + '" I= ' + str(N_i) + ' J= ' + str(
                    N_j) + ' K= ' + str(N_k) + '\n'

                # Create the block of data we want to write:
                data = np.zeros([N_i * N_j * N_k, N_var + 3])
                variables = self.rows[i].passage_instances[j].flow_field.variables

                data[:, 0] = np.reshape(self.rows[i].passage_instances[j].grid.x, -1, order='F')  # x
                data[:, 1] = np.reshape(self.rows[i].passage_instances[j].grid.y, -1, order='F')  # y
                data[:, 2] = np.reshape(self.rows[i].passage_instances[j].grid.z, -1, order='F')  # z

                for k in range(N_var):
                    data[:, 3 + k] = np.reshape(variables[computed_variables_names[k]], -1, order='F')

                # Write down the chink of data:
                np.savetxt(fname=file, X=data, fmt='%10.5e', header=header, comments='')

        file.close()
        print('Done! your file is "ParaView_TecPlotInterpreter_passages.dat"')
        print(' ')
