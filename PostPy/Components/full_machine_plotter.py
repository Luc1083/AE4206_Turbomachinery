import matplotlib.pyplot as plt
import numpy as np


class PlotMachine:
    def __init__(self, machine):
        self.machine = machine

    def convergence_history(self,file='stage.log'):
        try:
            file = open(file,'rt')
        except:
            print('The file specified or "stage.log" is not found, convergence history cannot be retrieved.')
            return
        data_residual = [line[15:68].split() for line in file]
        [data_residual[i].pop(3) for i in range(len(data_residual))]
        file.seek(0)
        data_pos = [line[79:].split() for line in file]
        file.close()

        N_data = len(data_pos)
        for i in range(N_data):
            data_residual[i][:] = [float(data_residual[i][j]) for j in range(4)]
            data_pos[i][:] = [int(data_pos[i][j]) for j in range(3)]
        data_residual = np.array(data_residual)
        data_pos = np.array(data_pos)
        data_pos = np.array(data_pos)
        iter = [i*5 for i in range(1,N_data+1)]

        N_grid = [self.machine.raw_data_grid[1][0], self.machine.raw_data_grid[1][1], self.machine.raw_data_grid[1][2]]

        fig, ax = plt.subplots(3)
        names = [r'$max\left| \frac{\partial u_m}{\partial t} \right|$', r'$avr\left| \frac{\partial u_m}{\partial t} \right|$', r'$max\left|\frac{\dot{m}}{\dot{m}_{in}}-1\right|$', 'FLOWRAT']
        for i in range(3):
            if i == 0:
                ax[0].plot(iter,[abs(data_residual[j,i]) for j in range(N_data)],label=names[i])
            else:
                ax[0].plot(iter, list(data_residual[:, i]), label=names[i])
        ax[0].set_yscale('log')
        ax[0].set_ylabel('Res')
        ax[0].legend()
        ax[0].set_title('Residuals')
        ax[0].set_xlim(left=0,right=(N_data+1)*5)

        names = [r'$\frac{I}{N_I}$', r'$\frac{J}{N_J}$', r'$\frac{K}{N_K}$']
        for i in range(3):
            ax[1].plot(iter, list(data_pos[:, i]/N_grid[i]), label=names[i])
        ax[1].set_title('Position of maximum error')
        ax[1].set_ylabel(r'index/$N_{max}$')
        ax[1].legend()
        ax[1].set_xlim(left=0,right=(N_data+1)*5)
        ax[1].set_ylim(0,1)

        ax[2].plot(iter,list(data_residual[:,3]))
        ax[2].set_title('Inlet mass flow')
        ax[2].set_ylabel('[Kg/s]')
        ax[2].set_xlim(left=0,right=(N_data+1)*5)
        plt.show()

    def variable_evolution_1D(self, variable, avr_type='massFlowAve'):
        N_rows = self.machine.N_rows
        N_k = self.machine.raw_data_grid[1][2]
        k = int(np.round(0.5*(N_k-1)))

        plt.figure()
        for r in range(N_rows):
            N_j = self.machine.rows[r].passage_original.N_j
            if avr_type == 'areaAve':
                avr = [self.machine.rows[r].passage_original.get_area_average(variable, j/(N_j-1), [0,1], [0,1]) for j in range(N_j)]
            elif avr_type == 'massFlowAve':
                avr = [self.machine.rows[r].passage_original.get_mass_flow_average(variable, j/(N_j-1), [0,1], [0,1]) for j in range(N_j)]
            else:
                raise TypeError("avr_type not recognised. it must be 'areaAve' or 'massFlowAve'")

            x = self.machine.rows[r].passage_original.grid.x[0,:,k].tolist()
            plt.plot(x,avr)

        plt.xlabel('x [m]')
        plt.ylabel(variable + ' [SI units]')
        plt.title(avr_type)
        plt.show()

    def variable_evolution_2D(self, variable, avr_type='massFlowAve'):
        N_rows = self.machine.N_rows

        plt.figure()
        var = []
        x = []
        y = []
        for r in range(N_rows):
            if (not avr_type == 'areaAve') and (not avr_type == 'massFlowAve'):
                raise TypeError("avr_type not recognised. it must be 'areaAve' or 'massFlowAve'")

            N_j = self.machine.rows[r].passage_original.N_j
            for j in range(N_j):
                var.append( self.machine.rows[r].passage_original.get_pitch_average(variable, j / (N_j-1), avr_type).tolist() )
                x.append(self.machine.rows[r].passage_original.grid.x[0,j,:].tolist())
                y.append(self.machine.rows[r].passage_original.grid.r[0, j, :].tolist())

            self.machine.rows[r].passage_original.grid.plot_contour_process('i', 0)
            self.machine.rows[r].blade_original.grid.plot_contour_process('i', 0,(1, 0, 0))

        var = np.array(var).transpose()
        x = np.array(x).transpose()
        y = np.array(y).transpose()

        plt.contourf(x,y,var)
        plt.colorbar()
        plt.xlabel('x [m]')
        plt.ylabel('r [m]')
        plt.suptitle(variable + ' [SI units], Pitch-wise ' + avr_type)
        plt.title('Black: row domains;  Red: blades contour')
        plt.axis('equal')
        plt.show()

    def variable_B2B(self,variable,k,levels):
        N_rows = self.machine.N_rows
        plt.figure()
        for i in range(N_rows):

            N_instances = self.machine.rows[i].N_instances
            if not N_instances == len(self.machine.rows[i].passage_instances):
                print('Updating geometry in row ' + str(i) + '...')
                self.machine.rows[i].generate_extra_geometry()
            for j in range(N_instances):
                self.machine.rows[i].passage_instances[j].plot_B2B_process(variable,k,levels)
                self.machine.rows[i].blade_instances[j].grid.plot_contour_process('k', k)
            self.machine.rows[i].blade_instances[N_instances].grid.plot_contour_process('k', k)

        plt.axis('equal')
        plt.colorbar()
        plt.xlabel('x [m]')
        plt.ylabel(r'$r\cdot\theta$ [m]')
        plt.title('Distribution of ' + variable + ' [SI units] at k=' + str(k) + '/ Span')

    def blades_contour(self,normal,level):
        N_rows = self.machine.N_rows
        plt.figure()
        for i in range(N_rows):
            self.machine.rows[i].blade_original.grid.plot_contour_process(normal, level)

    def blades_grid(self,normal,level):
        N_rows = self.machine.N_rows
        plt.figure()
        for i in range(N_rows):
            self.machine.rows[i].blade_original.grid.plot_grid_process(normal,level)

    def passage_contour(self, normal, level):
        N_rows = self.machine.N_rows
        plt.figure()
        for i in range(N_rows):
            self.machine.rows[i].passage_original.grid.plot_contour_process(normal, level)

    def passage_grid(self, normal, level):
        N_rows = self.machine.N_rows
        plt.figure()
        for i in range(N_rows):
            self.machine.rows[i].passage_original.grid.plot_grid_process(normal, level)

    def linear_cascade_blades(self, level):
        N_rows = self.machine.N_rows
        plt.figure()
        for i in range(N_rows):
            if not self.machine.rows[i].N_instances == len(self.machine.rows[i].passage_instances):
                print('Updating geometry in row ' + str(i) + '...')
                self.machine.rows[i].generate_extra_geometry()
            N_instances = self.machine.rows[i].N_instances+1
            for j in range(N_instances):
                self.machine.rows[i].blade_instances[j].grid.plot_contour_process('k',level)

    def linear_cascade_contour(self,level):
        N_rows = self.machine.N_rows
        plt.figure()
        for i in range(N_rows):
            if not self.machine.rows[i].N_instances == len(self.machine.rows[i].passage_instances):
                print('Updating geometry in row ' + str(i) + '...')
                self.machine.rows[i].generate_extra_geometry()
            N_instances = self.machine.rows[i].N_instances
            for j in range(N_instances):
                self.machine.rows[i].passage_instances[j].grid.plot_contour_process('k',level)

    def linear_cascade_grid(self,level):
        N_rows = self.machine.N_rows
        plt.figure()
        for i in range(N_rows):
            if not self.machine.rows[i].N_instances == len(self.machine.rows[i].passage_instances):
                print('Updating geometry in row ' + str(i) + '...')
                self.machine.rows[i].generate_extra_geometry()
            N_instances = self.machine.rows[i].N_instances
            for j in range(N_instances):
                self.machine.rows[i].passage_instances[j].grid.plot_grid_process('k',level)
