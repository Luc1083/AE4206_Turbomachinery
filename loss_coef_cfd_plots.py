import re
import numpy as np
import matplotlib.pyplot as plt
with open("results.out") as file:
    raw = file.readlines()
enthalpy_string=str(raw[-25])
temp_enthalpy = re.findall(r'[\d.]+', enthalpy_string)
enthalpy = list(map(float, temp_enthalpy))

massflow_string=str(raw[-26])
temp_massflow = re.findall(r'[\d.]+', massflow_string)
massflow_in, massflow_out = list(map(float, temp_massflow))

stag_temperature_string=str(raw[-27])
temp_stag_temperature = re.findall(r'[\d.]+', stag_temperature_string)
stag_temperature_in, stag_temperature_out  = list(map(float, temp_stag_temperature))

static_pres_string=str(raw[-28])
temp_static_pres = re.findall(r'[\d.]+', static_pres_string)
static_pres_in, static_pres_out = list(map(float, temp_static_pres))

stag_pres_string=str(raw[-29])
temp_stag_pres = re.findall(r'[\d.]+', stag_pres_string)
stag_pres_in, stag_pres_out = list(map(float, temp_stag_pres))

constants_string=str(raw[-35])
temp_constants = re.findall(r'[\d.]+', constants_string)
Cp, gamma, R = list(map(float, temp_constants))

temp=str(raw[-1])
temp_values = re.findall(r'[\d.]+', temp)
step,emax,I,J,K,eavg,econt,j2,vref,vmax,i2,j3,k3,mass_in,mass_out,rat_flow = list(map(float, temp_values))

Entropy_change=Cp*np.log(stag_temperature_out/stag_temperature_in)-R*np.log(stag_pres_out/stag_pres_in)
shock_loss=Entropy_change/(R*((static_pres_out-static_pres_in)/static_pres_in))

rho=0.55
print(vref)
BL_loss=(stag_pres_out-stag_pres_in)/(0.5*rho*vref**2)

loss_coef=[]
for i in range(-274,-52):
    loss_str=str(raw[i])
    loss_values = re.findall(r'[\d.]+', loss_str)
    i=+1
    loss=float(loss_values[3])
    loss_coef.append(loss)

domainlength=0.230919764+0.041058
domain=np.linspace(0,domainlength,len(loss_coef))


plt.figure()
plt.plot(domain,loss_coef)
plt.vlines(0.08012+0.041058,ymax=max(loss_coef),ymin=0,linestyles='dashed',colors='r')
plt.vlines(0.041058,ymax=max(loss_coef),ymin=0,linestyles='dashed',colors='r')
plt.vlines(0.19251+0.041058,ymax=max(loss_coef),ymin=0,linestyles='dashed',colors='y')
plt.vlines(0.11572+0.041058,ymax=max(loss_coef),ymin=0,linestyles='dashed',colors='y')
plt.xlabel('x [m]')
plt.ylabel('Lost efficiency')
plt.show()