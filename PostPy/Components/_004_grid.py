import numpy as np
import matplotlib.pyplot as plt


class Coordinates:
    def __init__(self,m,r,rth):
        self.m = m
        self.r = r
        self.rth = rth

        self.x = self.m
        self.y = self.r * np.sin(self.rth / self.r)
        self.z = self.r * np.cos(self.rth / self.r)

    def plot_contour(self,normal,level):
        plt.figure()
        self.plot_contour_process(normal,level)

    def plot_grid(self,normal,level):
        plt.figure()
        self.plot_grid_process(normal,level)

    def plot_contour_process(self,normal,level,colour=(0,0,0)):
        # Plots the outer bound of the grid in a plane perpendicular to 'normal' (i,I,j,J,k,K) at a level normal = l [0,1]
        if normal == 'i' or normal == 'I':
            normal_ind = 0
            plot_ind = [1,2]
            x = self.m
            y = self.r
            xlabel = 'Meridional direction [m]'
            ylabel = 'Radial direction [m]'
        elif normal == 'j' or normal == 'J':
            normal_ind = 1
            plot_ind = [0,2]
            x = self.z
            y = self.y
            xlabel = 'z [m]'
            ylabel = 'y [m]'
        elif normal == 'k' or normal == 'K':
            normal_ind = 2
            plot_ind = [1,0]
            x = self.m
            y = self.rth
            xlabel = 'Meridional direction [m]'
            ylabel = 'Tangential direction r·θ [m]'
        else:
            raise TypeError("The direction specified is not allowed, use one of these: i, I, j, J, k, K.")

        if level > 1 or level < 0:
            raise Exception("Level value must be between 0 and 1")

        size = self.m.shape
        N_n = size[normal_ind]
        N_x = size[plot_ind[0]]
        N_y = size[plot_ind[1]]

        l = round((N_n-1)*level)
        cords = [0,0,0]
        cords[normal_ind] = l

        cords[plot_ind[0]] = list(range(N_x))
        cords[plot_ind[1]] = 0
        cords = tuple(cords)
        line_1 = [x[cords],y[cords]]

        cords = list(cords)
        cords[plot_ind[0]] = list(range(N_x))
        cords[plot_ind[1]] = N_y-1
        cords = tuple(cords)
        line_2 = [x[cords], y[cords]]

        cords = list(cords)
        cords[plot_ind[0]] = 0
        cords[plot_ind[1]] = list(range(N_y))
        cords = tuple(cords)
        line_3 = [x[cords], y[cords]]

        cords = list(cords)
        cords[plot_ind[0]] = N_x-1
        cords[plot_ind[1]] = list(range(N_y))
        cords = tuple(cords)
        line_4 = [x[cords], y[cords]]

        # plt.figure()

        plt.plot(line_1[0], line_1[1], color=colour)
        plt.plot(line_2[0], line_2[1], color=colour)
        plt.plot(line_3[0], line_3[1], color=colour)
        plt.plot(line_4[0], line_4[1], color=colour)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.suptitle('Domain countour seen from direction ' + normal + '.')
        plt.title('% of displacement in that direction: ' + str(level))
        plt.axis('equal')

        plt.show()

    def plot_grid_process(self,normal,level):
        # Plots the outer bound of the grid in a plane perpendicular to 'normal' (i,I,j,J,k,K) at a level normal = l [0,1]
        if normal == 'i' or normal == 'I':
            normal_ind = 0
            plot_ind = [1,2]
            x = self.m
            y = self.r
            xlabel = 'Meridional direction [m]'
            ylabel = 'Radial direction [m]'
        elif normal == 'j' or normal == 'J':
            normal_ind = 1
            plot_ind = [0,2]
            x = self.z
            y = self.y
            xlabel = 'z [m]'
            ylabel = 'y [m]'
        elif normal == 'k' or normal == 'K':
            normal_ind = 2
            plot_ind = [1,0]
            x = self.m
            y = self.rth
            xlabel = 'Meridional direction [m]'
            ylabel = 'Tangential direction r·θ [m]'
        else:
            raise TypeError("The direction specified is not allowed, use one of these: i, I, j, J, k, K.")

        if level > 1 or level < 0:
            raise Exception("Level value must be between 0 and 1")

        size = self.m.shape
        N_n = size[normal_ind]
        N_x = size[plot_ind[0]]
        N_y = size[plot_ind[1]]

        l = round((N_n-1)*level)
        cords = [0,0,0]
        cords[normal_ind] = l

        # plt.figure()

        for x_ in range(N_x):
            cords = list(cords)
            cords[plot_ind[0]] = x_
            cords[plot_ind[1]] = list(range(N_y))
            cords = tuple(cords)
            line_1 = [x[cords],y[cords]]
            plt.plot(line_1[0], line_1[1], 'k')

        for y_ in range(N_y):
            cords = list(cords)
            cords[plot_ind[0]] = list(range(N_x))
            cords[plot_ind[1]] = y_
            cords = tuple(cords)
            line_1 = [x[cords],y[cords]]
            plt.plot(line_1[0], line_1[1], 'k')

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.suptitle('Grid seen from direction ' + normal + ' at %: ' + str(level) + '.')
        plt.axis('equal')

        plt.show()


