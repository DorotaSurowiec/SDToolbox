# Updated 2-20-2014 for compitability with
# Python 2.7 and Cantera 2.1.
import sys
import numpy as np
from cantera import *
from matplotlib.pylab import *
import csv
from SDToolbox import *

#Keybord input
Tmin = 1000
Tmax = 2000
Pmin = 101325
Pmax = 405300
fimin = 0.2
fimax = 3.8

#Number of iterations
npoints = 11
fipoints = 10

#Creating lists for data storage
Pi = np.zeros(npoints, 'd')
Ti = np.zeros(npoints, 'd')
fi = np.zeros(fipoints, 'd')

#Lists for plots
cj_speed_casPressure = np.zeros(npoints, 'd')
rho_casPressure = np.zeros(npoints, 'd')
press_casPressure = np.zeros(npoints, 'd')
temp_casPressure = np.zeros(npoints, 'd')
af_casPressure = np.zeros(npoints, 'd')
ae_casPressure = np.zeros(npoints, 'd')

cj_speed_casTemp = np.zeros(npoints, 'd')
rho_casTemp = np.zeros(npoints, 'd')
press_casTemp = np.zeros(npoints, 'd')
temp_casTemp = np.zeros(npoints, 'd')
af_casTemp = np.zeros(npoints, 'd')
ae_casTemp = np.zeros(npoints, 'd')

cj_speed_casFi = np.zeros(fipoints, 'd')
rho_casFi = np.zeros(fipoints, 'd')
press_casFi = np.zeros(fipoints, 'd')
temp_casFi = np.zeros(fipoints, 'd')
af_casFi = np.zeros(fipoints, 'd')
ae_casFi = np.zeros(fipoints, 'd')

#Ethane-Oxygen
#T=const and phi=const
for p in range(npoints):
    Pi[p] = Pmin + (Pmax-Pmin)*p/(npoints-1)
    T1=1000;
    fi1=1
    q='C2H6:0.2857 O2:1';
    mech = 'gri30_highT.cti'; #mechanism used for the process
 
    [cj_speed,R2] = CJspeed(Pi[p], T1, q, mech, 0);   

    gas = PostShock_eq(cj_speed, Pi[p], T1, q, mech)
    Ps = gas.P/one_atm

    gas()
    
    [ae,af] = equilSoundSpeeds(gas)
    
    cj_speed_casPressure[p] = cj_speed
    rho_casPressure[p] = gas.density
    press_casPressure[p] = gas.P/one_atm
    temp_casPressure[p] = gas.T
    af_casPressure[p] = af
    ae_casPressure[p] = ae
    
csv_file='sdt_casPressure.csv'
with open(csv_file,'w') as outfile:
    writer=csv.writer(outfile)
    writer.writerow(['Initial temperature','Pressure','phi','cjspeed','Final Temperature', 'Final pressure', 'density', 'af speed', 'ae speed'])
    for i in range(npoints):
        writer.writerow([T1,Pi[i],fi1,cj_speed_casPressure[i],temp_casPressure[i],press_casPressure[i],rho_casPressure[i], af_casPressure[i], ae_casPressure[i]])
        
#P=const and phi=const
for t in range(npoints):
    Ti[t] = Tmin + (Tmax-Tmin)*t/(npoints-1)
    P1=101325;
    fi1=1
    q='C2H6:0.2857 O2:1';
    mech = 'gri30_highT.cti'; #mechanism used for the process

 
    [cj_speed,R2] = CJspeed(P1, Ti[t], q, mech, 0);   

    gas = PostShock_eq(cj_speed, P1, Ti[t], q, mech)
    Ps = gas.P/one_atm

    gas()
    
    [ae,af] = equilSoundSpeeds(gas)
    
    cj_speed_casTemp[t] = cj_speed
    rho_casTemp[t] = gas.density
    press_casTemp[t] = gas.P/one_atm
    temp_casTemp[t] = gas.T
    af_casTemp[t] = af
    ae_casTemp[t] = ae

csv_file='sdt_casTemp.csv'
with open(csv_file,'w') as outfile:
    writer=csv.writer(outfile)
    writer.writerow(['Initial temperature','Pressure','phi','cjspeed','Final Temperature', 'Final pressure', 'density', 'af speed', 'ae speed'])
    for i in range(npoints):
        writer.writerow([Ti[i],P1,fi1,cj_speed_casTemp[i],temp_casTemp[i],press_casTemp[i],rho_casTemp[i], af_casTemp[i], ae_casTemp[i]])

#T=const and P=const
for f in range(fipoints):
    fi[f] = fimin + (fimax-fimin)*f/(fipoints-1)
    P1=101325;
    T1=1000
    no = float(1/fi[f])
    q='C2H6:0.2857 O2:' + str(no);
    mech = 'gri30_highT.cti'; #mechanism used for the process

 
    [cj_speed,R2] = CJspeed(P1, T1, q, mech, 0);   

    gas = PostShock_eq(cj_speed, P1, T1, q, mech)
    Ps = gas.P/one_atm

    gas()
    
    [ae,af] = equilSoundSpeeds(gas)
    
    cj_speed_casFi[f] = cj_speed
    rho_casFi[f] = gas.density
    press_casFi[f] = gas.P/one_atm
    temp_casFi[f] = gas.T
    af_casFi[f] = af
    ae_casFi[f] = ae

csv_file='sdt_casFi.csv'
with open(csv_file,'w') as outfile:
    writer=csv.writer(outfile)
    writer.writerow(['Initial temperature','Pressure','phi','cjspeed','Final Temperature', 'Final pressure', 'density', 'af speed', 'ae speed'])
    for i in range(fipoints):
        writer.writerow([T1,P1,fi[i],cj_speed_casFi[i],temp_casFi[i],press_casFi[i],rho_casFi[i], af_casFi[i], ae_casFi[i]])

#PlOTS
#for pressure
plot(Pi,cj_speed_casPressure,'-',color='blue')
xlabel(r'Pressure [Pa]',fontsize=20)
ylabel("CJ speed [m/s]")
title(r'CJ speed', fontsize=22,horizontalalignment='center')
axis([100000,410000,2200,2400])
grid()
savefig('press_cjspeed.png',bbox_inches='tight')

plot(Pi,rho_casPressure,'-',color='blue')
xlabel(r'Pressure [Pa]',fontsize=20)
ylabel("Density [kg/m3]")
title(r'Density', fontsize=22,horizontalalignment='center')
axis([100000,410000,0.5,3.0])
grid()
savefig('press_rho.png',bbox_inches='tight')

plot(Pi,press_casPressure,'-',color='blue')
xlabel(r'Pressure [Pa]',fontsize=20)
ylabel("Pressure [Pa]")
title(r'Final pressure', fontsize=22,horizontalalignment='center')
axis([100000,410000,5,45])
grid()
savefig('press_pressure.png',bbox_inches='tight')

plot(Pi,temp_casPressure,'-',color='blue')
xlabel(r'Pressure [Pa]',fontsize=20)
ylabel("Temperature [K]")
title(r'Final Temperature', fontsize=22,horizontalalignment='center')
axis([100000,410000,3600,4000])
grid()
savefig('press_temp.png',bbox_inches='tight')

plot(Pi,af_casPressure,'b-',Pi,ae_casPressure,'r-')
xlabel(r'Pressure [Pa]',fontsize=20)
ylabel("speed of sound [m/s]")
title(r'Af and ae speed of sound', fontsize=22,horizontalalignment='center')
axis([100000,410000,1280,1380])
grid()
savefig('press_soundspeed.png',bbox_inches='tight')

#for temperature

plot(Ti,cj_speed_casTemp,'-',color='blue')
xlabel(r'Temperature [K]',fontsize=20)
ylabel("CJ speed [m/s]")
title(r'CJ speed', fontsize=22,horizontalalignment='center')
axis([1000,2000,2200,2350])
grid()
savefig('temp_cjspeed.png',bbox_inches='tight')

plot(Ti,rho_casTemp,'-',color='blue')
xlabel(r'Temperature [K]',fontsize=20)
ylabel("Density [kg/m3]")
title(r'Density', fontsize=22,horizontalalignment='center')
axis([1000,2000,0,1])
grid()
savefig('temp_rho.png',bbox_inches='tight')

plot(Ti,press_casTemp,'-',color='blue')
xlabel(r'Temperature [K]',fontsize=20)
ylabel("Pressure [Pa]")
title(r'Final pressure', fontsize=22,horizontalalignment='center')
axis([1000,2000,5,10])
grid()
savefig('temp_pressure.png',bbox_inches='tight')

plot(Ti,temp_casTemp,'-',color='blue')
xlabel(r'Temperature [K]',fontsize=20)
ylabel("Temperature [K]")
title(r'Final Temperature', fontsize=22,horizontalalignment='center')
axis([1000,2000,3600,3700])
grid()
savefig('temp_temp.png',bbox_inches='tight')

plot(Ti,af_casTemp,'b-',Ti,ae_casTemp,'r-')
xlabel(r'Temperature [K]',fontsize=20)
ylabel("speed of sound [m/s]")
title(r'Af and ae speed of sound', fontsize=22,horizontalalignment='center')
axis([1000,2000,1280,1410])
grid()
savefig('temp_soundspeed.png',bbox_inches='tight')

#for phi
plot(fi,cj_speed_casFi,'-',color='blue')
xlabel(r'phi',fontsize=20)
ylabel("CJ speed [m/s]")
title(r'CJ speed', fontsize=22,horizontalalignment='center')
axis([0.0,4.0,1650,2650])
grid()
savefig('phi_cjspeed.png',bbox_inches='tight')

plot(fi,rho_casFi,'-',color='blue')
xlabel(r'phi',fontsize=20)
ylabel("Density [kg/m3]")
title(r'Density', fontsize=22,horizontalalignment='center')
axis([0.0,4.0,0.0,1.0])
grid()
savefig('phi_rho.png',bbox_inches='tight')


plot(fi,press_casFi,'-',color='blue')
xlabel(r'phi',fontsize=20)
ylabel("Pressure [Pa]")
title(r'Final pressure', fontsize=22,horizontalalignment='center')
axis([0.0,4.0,5,15])
grid()
savefig('phi_pressure.png',bbox_inches='tight')

plot(fi,temp_casFi,'-',color='blue')
xlabel(r'phi',fontsize=20)
ylabel("Temperature [K]")
title(r'Final Temperature', fontsize=22,horizontalalignment='center')
axis([0.0,4.0,2000,3700])
grid()
savefig('phi_temp.png',bbox_inches='tight')

plot(fi,af_casFi,'b-',fi,ae_casFi,'r-')
xlabel(r'phi',fontsize=20)
ylabel("speed of sound [m/s]")
title(r'Af and ae speed of sound', fontsize=22,horizontalalignment='center')
axis([0.0,4.0,950,1550])
grid()
savefig('phi_soundspeed.png',bbox_inches='tight')