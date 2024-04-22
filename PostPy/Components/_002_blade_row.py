import numpy as np
from ._003_multiblock_solution import BlockCFD
from copy import deepcopy


class BladeRow:
    def __init__(self, grid_data, flow_data,ref_values):
        self.raw_data_grid = grid_data
        self.raw_data_flow = flow_data
        self.reference_values = ref_values

        self.J_TE = None
        self.J_LE = None

        if self.raw_data_grid[4][1] == 0:
            self.type = 'stator'
        else:
            self.type = 'rotor'

        self.N_blades = self.raw_data_grid[5][1]
        self.RPM = self.raw_data_grid[4][1] * 60 / (2 * np.pi)
        [data_blade_grid, data_blade_flow] = self.extract_blade()

        self.passage_original = BlockCFD(self.raw_data_grid, self.raw_data_flow, ref_values)  # Original copy of the passage
        self.blade_original = BlockCFD(data_blade_grid, data_blade_flow, ref_values)  # Original copy of the blade

        self.N_instances = 2  # Number of copies of the passage that are going to be generated (default value)
        self.passage_instances = []  # Displaced copies of the passage itself (list)
        self.blade_instances = []  # Displaced copies of the blades (list)

    def extract_blade(self):
        # Size of the grid:
        N_i = self.raw_data_grid[1][0]
        N_j = self.raw_data_grid[1][1]
        N_k = self.raw_data_grid[1][2]

        # Extract blades
        # Find LE and TE
        flag = False
        for j in range(N_j):
            if self.raw_data_grid[3][j] == 1 and not flag:
                flag = True
                J_LE = j
            elif self.raw_data_grid[3][j] == 1 and flag:
                J_TE = j
        N_j_blade = J_TE-J_LE+1

        # Geometry (grid)
        delta_theta = 2 * np.pi / self.raw_data_grid[5][0]
        coordinates = self.raw_data_grid[6 + N_k*J_LE : 6 + N_k*(J_TE+1)]
        coordinates = [[item[0],item[1],item[2],item[N_i+1]] for item in coordinates]

        data_blade_grid = []
        data_blade_grid.append(self.raw_data_grid[0])
        data_blade_grid.append([2, N_j_blade, self.raw_data_grid[1][2]])
        data_blade_grid.append(self.raw_data_grid[2])
        for index in range(3, 6):
            data_blade_grid.append(self.raw_data_grid[index][J_LE:J_TE + 1])
        for i in range(len(coordinates)):
            coordinates[i][3] -= delta_theta * coordinates[i][1]
            data_blade_grid.append(coordinates[i])

        # Flow field:
        data_blade_flow = []
        data_blade_flow.append(self.raw_data_flow[0])
        for i_var in range(1,len(self.raw_data_flow)):
            aux_property = []
            for k in range(N_k):
                for j in range(J_LE,J_TE+1):
                    # The data is stored in a vector such that the i,k,k correspondence is index=i+j·Ni+k·(Ni·Nj), with the indices starting at 0 and finishing at N-1
                    aux_property.append(self.raw_data_flow[i_var][0 + j*N_i + k*(N_i*N_j)])
                    aux_property.append(self.raw_data_flow[i_var][(N_i-1) + j*N_i + k*(N_i*N_j)])
            data_blade_flow.append(aux_property)

        # Store values in class attributes:
        self.J_LE = J_LE
        self.J_TE = J_TE

        return data_blade_grid, data_blade_flow

    def generate_extra_geometry(self):
        self.blade_instances = []
        self.passage_instances = []
        delta_theta = 2*np.pi/self.N_blades
        for i in range(self.N_instances):
            # Create more passages:
            data_grid = deepcopy(self.passage_original.raw_data_grid)
            data_flow = deepcopy(self.passage_original.raw_data_flow)

            for j in range(6, len(data_grid)):
                for k in range(2, len(data_grid[j])):
                    data_grid[j][k] += i * delta_theta * data_grid[j][1]
            self.passage_instances.append(BlockCFD(data_grid,data_flow,self.reference_values))

            # Create more blades:
            data_grid = deepcopy(self.blade_original.raw_data_grid)
            data_flow = deepcopy(self.blade_original.raw_data_flow)
            for j in range(6, len(data_grid)):
                for k in range(2, len(data_grid[j])):
                    data_grid[j][k] += i * delta_theta * data_grid[j][1]
            self.blade_instances.append(BlockCFD(data_grid,data_flow,self.reference_values))

        # This is one extra blade outside the loop (N passages delimited by N+1 blades):
        data_grid = deepcopy(self.blade_original.raw_data_grid)
        for j in range(6, len(data_grid)):
            for k in range(2, len(data_grid[j])):
                data_grid[j][k] += self.N_instances * delta_theta * data_grid[j][1]

        self.blade_instances.append(BlockCFD(data_grid,data_flow,self.reference_values))
