import PostPy.Components as PP
import numpy as np

tm = PP.Turbomachine()
#tm.gen_ParaView_input()
#tm.plot.convergence_history()
#tm.rows[0].N_instances = 5
#tm.rows[1].N_instances = 5
#tm.plot.variable_B2B("M_rel", 0.99, 10)
#tm.plot.variable_B2B("P", 0.9, 10)
#tm.rows[0].blade_original.plot_Cp(0.5)
#tm.rows[0].blade_original.plot_pitch_average_evolution('s','massFlowAve')
#tm.plot.variable_evolution_2D('s') 
U_mid = tm.rows[0].passage_original.get_mass_flow_average('U', 0.5, [0,1], [0.95,1])
V_ax = tm.rows[0].passage_original.get_mass_flow_average('Vm_rel', 0.5, [0,1], [0,1])
w=3130.71704E3/77.1044617
flow=V_ax/U_mid
work=w/U_mid**2

rotor_loss=0.132505-0.0146701
stator_loss=0.184896-0.163478


s_tip_1 = tm.rows[0].passage_original.get_mass_flow_average('s', 0, [0,1], [0.95,1])
s_tip_2 = tm.rows[0].passage_original.get_mass_flow_average('s', 1, [0,1], [0.95,1])
s_tip_3 = tm.rows[1].passage_original.get_mass_flow_average('s', 1, [0,1], [0.95,1])
T_2_tip = tm.rows[0].passage_original.get_mass_flow_average('T', 1, [0,1], [0.95,1])
T_3_tip= tm.rows[1].passage_original.get_mass_flow_average('T', 1, [0,1], [0.95,1])
vx_tip=tm.rows[0].passage_original.get_mass_flow_average('Vx_stn', 1, [0,1], [0.95,1])
vy_tip=tm.rows[0].passage_original.get_mass_flow_average('Vy_stn', 1, [0,1], [0.95,1])
vz_tip=tm.rows[0].passage_original.get_mass_flow_average('Vz_stn', 1, [0,1], [0.95,1])
v2_tip=np.sqrt(vx_tip**2+vy_tip**2+vz_tip**2)
wx_tip=tm.rows[0].passage_original.get_mass_flow_average('Vx_rel', 1, [0,1], [0.95,1])
wy_tip=tm.rows[0].passage_original.get_mass_flow_average('Vy_rel', 1, [0,1], [0.95,1])
wz_tip=tm.rows[0].passage_original.get_mass_flow_average('Vz_rel', 1, [0,1], [0.95,1])
w1_tip=np.sqrt(wx_tip**2+wy_tip**2+wz_tip**2)

eta_rotor_tip= (T_2_tip*(abs(s_tip_2-s_tip_1)))/(0.5*w1_tip**2)
eta_stator_tip= (T_3_tip*(abs(s_tip_3-s_tip_2)))/(0.5*v2_tip**2)

s_hub_1 = tm.rows[0].passage_original.get_mass_flow_average('s', 0, [0,1], [0,0.05])
s_hub_2 = tm.rows[0].passage_original.get_mass_flow_average('s', 1, [0,1], [0,0.05])
s_hub_3 = tm.rows[1].passage_original.get_mass_flow_average('s', 1, [0,1], [0,0.05])
T_2_hub = tm.rows[0].passage_original.get_mass_flow_average('T', 1, [0,1], [0,0.05])
T_3_hub= tm.rows[1].passage_original.get_mass_flow_average('T', 1, [0,1], [0,0.05])
vx_hub=tm.rows[0].passage_original.get_mass_flow_average('Vx_stn', 1, [0,1], [0,0.05])
vy_hub=tm.rows[0].passage_original.get_mass_flow_average('Vy_stn', 1, [0,1], [0,0.05])
vz_hub=tm.rows[0].passage_original.get_mass_flow_average('Vz_stn', 1, [0,1], [0,0.05])
v2_hub=np.sqrt(vx_hub**2+vy_hub**2+vz_hub**2)
wx_hub=tm.rows[0].passage_original.get_mass_flow_average('Vx_rel', 1, [0,1], [0,0.05])
wy_hub=tm.rows[0].passage_original.get_mass_flow_average('Vy_rel', 1, [0,1], [0,0.05])
wz_hub=tm.rows[0].passage_original.get_mass_flow_average('Vz_rel', 1, [0,1], [0,0.05])
w1_hub=np.sqrt(wx_hub**2+wy_hub**2+wz_hub**2)

eta_rotor_hub= (T_2_hub*(abs(s_hub_2-s_hub_1)))/(0.5*w1_hub**2)
eta_stator_hub= (T_3_hub*(abs(s_hub_3-s_hub_2)))/(0.5*v2_hub**2)

s_mid_1 = tm.rows[0].passage_original.get_mass_flow_average('s', 0, [0,1], [0.45,0.55])
s_mid_2 = tm.rows[0].passage_original.get_mass_flow_average('s', 1, [0,1], [0.45,0.55])
s_mid_3 = tm.rows[1].passage_original.get_mass_flow_average('s', 1, [0,1], [0.45,0.55])
T_2_mid = tm.rows[0].passage_original.get_mass_flow_average('T', 1, [0,1], [0.45,0.55])
T_3_mid= tm.rows[1].passage_original.get_mass_flow_average('T', 1, [0,1], [0.45,0.55])
vx_mid=tm.rows[0].passage_original.get_mass_flow_average('Vx_stn', 1, [0,1], [0.45,0.55])
vy_mid=tm.rows[0].passage_original.get_mass_flow_average('Vy_stn', 1, [0,1], [0.45,0.55])
vz_mid=tm.rows[0].passage_original.get_mass_flow_average('Vz_stn', 1, [0,1], [0.45,0.55])
v2_mid=np.sqrt(vx_mid**2+vy_mid**2+vz_mid**2)
wx_mid=tm.rows[0].passage_original.get_mass_flow_average('Vx_rel', 1, [0,1], [0.45,0.55])
wy_mid=tm.rows[0].passage_original.get_mass_flow_average('Vy_rel', 1, [0,1], [0.45,0.55])
wz_mid=tm.rows[0].passage_original.get_mass_flow_average('Vz_rel', 1, [0,1], [0.45,0.55])
w1_mid=np.sqrt(wx_mid**2+wy_mid**2+wz_mid**2)

eta_rotor_mid= (T_2_mid*(abs(s_mid_2-s_mid_1)))/(0.5*w1_mid**2)
eta_stator_mid= (T_3_mid*(abs(s_mid_3-s_mid_2)))/(0.5*v2_mid**2)
























print(w1_hub)
print(v2_hub)

print(eta_rotor_tip)
print(eta_stator_tip)

eta=1-(1/(2*0.6))*((eta_rotor_hub*w1_hub**2+eta_stator_hub*v2_hub**2)/U_mid**2)
print(eta_rotor_hub)
print(eta_stator_hub)

print(eta_rotor_mid)
print(eta_stator_mid)
print(eta)


















