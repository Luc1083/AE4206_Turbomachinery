from Components import *

machine = Turbomachine('grid_out', 'flow_out')
N_rows = machine.N_rows

factor = 0.25
for i in range(N_rows):
    machine.rows[i].N_instances = int(np.ceil(machine.rows[i].N_blades * factor))

machine.gen_ParaView_input()
