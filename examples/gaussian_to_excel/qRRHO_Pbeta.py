#!/usr/bin/env python

# The format to execute this file is qRRO_Pbeta.py #1 #2 #3. #1 is the filename, #2 and #3 are desired T & P.
import sys
import os
import ast
import shutil as sh
from argparse import ArgumentParser as AP
import re
import numpy as np
import math


parser = AP(description='Extract thermal correction from Gaussain output' )


parser.add_argument('filename', type=str, default=None,
                    help='Name of log file from Gaussian')

parser.add_argument('-t', '--temperature',
                    default='180', type=float,
                    help='Temperature for thermal correction')

parser.add_argument('-p', '--pressure',default=1.0,
                    type=float, help='Pressure for thermal correction')

#parser.add_argument('-a', '--adsorbate',default=None, type=str,nargs='+', help='Name of the adsorbed molecules')


parser.add_argument('-n', '--number',default=0, type=int, help='Number of adsorbed molecules on the zeolite surface')

parser.add_argument('-m1', '--mass1',default=0,
                    type=int, help='Mass of the first adsorbate in atomic unit')

parser.add_argument('-m2', '--mass2',default=0,
                    type=int, help='Total mass of the second adsorbate in atomic unit')

parser.add_argument('-f', '--frequency',default=1.0,
                    type=float, help='frequency scale factor')

args = parser.parse_args()

file = args.filename
f = open(file,'r')
F=[]
#E = []
for line in f:
        if re.search('Zero-point correction=(.*?)\(', line):
            z = float(re.search('Zero-point correction=(.*?)\(', line).groups()[0])

        if re.search('Sum of electronic and zero-point Energies=(.*)',line):
            epz = float(re.search('Sum of electronic and zero-point Energies=(.*)',line).groups()[0])
        
       # if re.search('Sum of electronic and thermal Energies=(.*)',line):
       #     ept = float(re.search('Sum of electronic and thermal Energies=(.*)',line).groups()[0])
       # 
       # if re.search('Sum of electronic and thermal Enthalpies=(.*)',line):
       #     eph = float(re.search('Sum of electronic and thermal Enthalpies=(.*)',line) .groups()[0])
       # 
       # if re.search('Sum of electronic and thermal Free Energies=(.*)',line):
       #     epg = float(re.search('Sum of electronic and thermal Free Energies=(.*)',line).groups()[0])

    # Get all the frequencies in the unit of cm^-1.
        if re.search('Frequencies -- (.*)',line):

            a = re.search('Frequencies -- (.*)',line)
            l = (a.groups()[0]).split()
            #l = list(map(float,l))
            l = [float(i) for i in l]
            #print l.split()
            for i in range(0,len(l)):
                     F.append(l[i])
    
        w_0 = 1
        #if re.search('Charge =  0 Multiplicity =(.*)',line):
        #    w_0 = float(re.search('Charge =  0 Multiplicity =(.*)',line) .groups()[0])
        
        if re.search('Rotational temperatures \(Kelvin\)(.*)',line):
            R_t = (re.search('Rotational temperatures \(Kelvin\)(.*)',line) .groups()[0]).split()   
            R_t = [float(i) for i in R_t]
            #print (R_t)

        if re.search('Molecular mass:(.*)',line):
            amu = (re.search('Molecular mass:(.*)',line) .groups()[0]).split()   
            amu = float(amu[0])
            #print (amu)

        if re.search('Rotational symmetry number(.*)',line):
            sigma = (re.search('Rotational symmetry number(.*)',line) .groups()[0]).split()   
            sigma = float(sigma[0])
            
                       
     # Get the individual contribution to the internal thermal energy at new temperature(E_tot = E_e+E_t+E_r+E_v)

#Get the desired temperature and pressure from input 
T = args.temperature + 273.15 #convert C to K
P = args.pressure * 101325 #Pa
c = 29979245800 #cm/s
pi = math.pi
k_B = 1.3806488 * 10**-23 #J/K
m = amu * 1.660539040* 10**-27 #kg
h = 6.626070040 * 10**-34 #j s
R = 8.3144598   #j/mol K

# Translational partition, energy, entropy
q_t = (2*pi*m)**1.5 * (k_B*T)**2.5 / (h**3 *P)
E_t = 1.5*R*T #J/mol
S_t = R*(math.log(q_t)+2.5) #J/mol K

# 2D translational motion for adsorbate


#def g(x):
#    return {
#    'AA' : 102 * 1.660539040* 10**-27 * 2, #kg mass of the adsorbate AA
#    'MF_1': 82 * 1.660539040* 10**-27 *2, #kg mass of the adsorbate MF
#    'MF' : 184 * 1.660539040* 10**-27, #kg mass of the adsorbate MF and AA
#    'DMF' : 198 * 1.660539040* 10**-27, #kg mass of the adsorbate DMF and AA
#    }[x]


#if args.mass1:
m_ads1 = args.mass1 * 1.660539040* 10**-27 #kg mass of the first adsorbat
m_ads2 = args.mass2 * 1.660539040* 10**-27 #kg mass of the first adsorbat    
n_ads = args.number #get the number of the adsorbed molecule
# 2 Dimension
A = 800**2 * 10**-24 #m^2 surface area of BEA
#A = 600 * 200 * 10**-24 #m^2 surface area of ZSM5
q_t_2d_1 = (2*pi*m_ads1*k_B*T)/(h**2) * A 
q_t_2d_2 = (2*pi*m_ads2*k_B*T)/(h**2) * A
if args.mass1 == 0 :
    q_t_2d_1 = 1
if args.mass2 == 0 :
    q_t_2d_2 = 1  
q_t_2d = q_t_2d_1 * q_t_2d_2 
E_t_2d = n_ads*R*T  #J/mole
S_t_2d = R*(math.log(q_t_2d)+n_ads)  #J/mol K
# 1 Dimendion
L = 800 * 10**-12 #m distance of BEA
#L = 600 * 10**-12 #m distance of ZSM5
q_t_1d_1 = math.sqrt(2*pi*m_ads1*k_B*T/h**2)*L
q_t_1d_2 = math.sqrt(2*pi*m_ads2*k_B*T/h**2)*L
if args.mass1 == 0 :
    q_t_1d_1 = 1
if args.mass2 == 0 :
    q_t_1d_2 = 1
q_t_1d = q_t_1d_1 * q_t_1d_2
E_t_1d = n_ads*0.5*R*T #J/mole
S_t_1d = R*(math.log(q_t_1d)+n_ads*0.5)   #J/mol K

# Electronic partition, energy, entropy
q_e = w_0
E_e = 0
S_e = R*(math.log(q_e))

# Rotational partition, energy, entropy
q_r = math.sqrt(pi/(R_t[0]*R_t[1]*R_t[2])) * T**1.5/sigma
E_r = 1.5*R*T
S_r = R*(math.log(q_r)+1.5)

#print (q_r)
# Vibrational partition, energy, entropy

# Check if there is any imaginary frequency
if F[0]<0 :
    n_i = 1 #indicator of one imaginary frequency
    F.remove(F[0])

f_scale = args.frequency #get the scale factor of the frequency

F = [i*f_scale for i in F]

# Vibrational temperature
theta = [h*i*c/k_B for i in F]
# frequency in s^-1
niu = [i*c for i in F]
#print (theta)

V = [R*i*(0.5+1/(math.exp(i/T)-1)) for i in theta]
E_v = R*sum(i*(0.5+1/(math.exp(i/T)-1)) for i in theta)
S_v_HO = [R*(i/(T*(math.exp(i/T)-1))-math.log(1-math.exp(-i/T))) for i in theta] #Consider vibration as harmonic oscillation
#S_v = R*sum(i/(T*(math.exp(i/T)-1))-math.log(1-math.exp(-i/T)) for i in theta)
S_v = sum(S_v_HO)

q_v = 1
for i in theta:
    q_v = q_v * math.exp(-i/(2*T))/(1-math.exp(-i/T)) 

#Using qRRHO approximation
w = [1/(1+(100/i)**4) for i in F]
E_v_qRRHO = sum(i*j + (1-i)*0.5*R*T for i,j in zip(w,V))

mu = [h/(8*pi**2*i) for i in niu] #moment of inertia for a free rotor with the same frequency
B_av = 1e-44 #kg m^2 Average molecular moment of inertia as a limiting value for small niu
mu_p = [i*B_av/(i+B_av) for i in mu]
#mu_p = mu
#print ([i - j for i, j in zip(mu,mu_p)])
S_v_r = [R*(0.5 + math.log(math.sqrt(8*pi**3*i*k_B*T/h**2))) for i in mu_p]
S_v_r_qRRHO = [i*j + (1-i)*k for i,j,k in zip(w,S_v_HO,S_v_r)]
S_v_qRRHO = sum (S_v_r_qRRHO)

#print (S_v/4184, S_v_qRRHO/4184)
# Total thermal correction
E_tot = (E_t + E_r + E_v + E_e)/4184 #kcal/mol

E_tot_qRRHO = (E_t + E_r + E_v_qRRHO + E_e)/4184 #kcal/mol
# Thermal correction to enthalpy

H_corr = E_tot + R*T/4184 #kcal/mol
H_corr_qRRHO = E_tot_qRRHO + R*T/4184 #kcal/mol

# Total Entropy
S_tot = (S_t + S_r +S_v +S_e)/4184 #kcal/mol K
S_tot_qRRHO = (S_t + S_r +S_v_qRRHO +S_e)/4184 #kcal/mol K


# Thermal correction to Gibbs Free Energy
G_corr = H_corr - T*S_tot #kcal/mol
G_corr_qRRHO = H_corr_qRRHO - T*S_tot_qRRHO #kcal/mol

#save all the correction
Harm_corr = [E_tot,H_corr,G_corr]  #Thermal, enthalpy and free energy correction for harmonic approximation with all degrees of freedom
q_RRHO_corr = [E_tot_qRRHO,H_corr_qRRHO,G_corr_qRRHO]


#Electronic energy
e = (epz - z)*627.509469 #kcal/mol

# Energies at new T and P

N_epz = epz*627.509469
N_ept = e + E_tot
N_eph = e + H_corr
N_epg = e + G_corr

# Energies at new T and P with qRRHO
N_ept_qRRHO = e + E_tot_qRRHO
N_eph_qRRHO = e + H_corr_qRRHO
N_epg_qRRHO = e + G_corr_qRRHO

N_p = [e,z*627.509469,N_epz,N_ept,N_eph,N_epg,S_tot]
N_p_qRRHO = [N_epz,N_ept_qRRHO,N_eph_qRRHO,N_epg_qRRHO,S_tot_qRRHO]

# For solid or adsorbed state, only vibrational degrees of freedom are considered

S_E_tot = E_v/4184 #S stands for Solid
S_H_corr = S_E_tot+ R*T/4184
S_S_tot = (S_v+R)/4184 #entropy due to vibration. R comes from the Sterling approximation.
S_G_corr = S_H_corr - T*S_S_tot

S_ept = e + S_E_tot
S_eph = e + S_H_corr 
S_epg = e + S_G_corr

S_p = [N_epz,S_ept,S_eph,S_epg]

#No q_RRHO with 2D translational motion
#if args.mass1:
E_tot_2d = (E_v + E_t_2d)/4184 
H_corr_2d = E_tot_2d+ R*T/4184
S_tot_2d = (S_v + S_t_2d + R)/4184 #eneropy due to vibration and 2d translation
G_corr_2d = H_corr_2d - T*S_tot_2d
ept_2d = e + E_tot_2d
eph_2d = e + H_corr_2d
epg_2d = e + G_corr_2d
S_p_2d = [N_epz,ept_2d,eph_2d,epg_2d,S_tot_2d]

#No q_RRHO with 1D translational motion
#if args.mass1:
E_tot_1d = (E_v + E_t_1d)/4184 
H_corr_1d = E_tot_1d+ R*T/4184
S_tot_1d = (S_v + S_t_1d + R)/4184 #eneropy due to vibration and 1d translation
G_corr_1d = H_corr_1d - T*S_tot_1d
ept_1d = e + E_tot_1d
eph_1d = e + H_corr_1d
epg_1d = e + G_corr_1d
S_p_1d = [N_epz,ept_1d,eph_1d,epg_1d,S_tot_1d]
 
#q_RRHO without any translational motion
S_E_tot_qRRHO = E_v_qRRHO/4184 
S_H_corr_qRRHO = S_E_tot_qRRHO+ R*T/4184
S_S_tot_qRRHO = (S_v_qRRHO + R)/4184 #eneropy due to vibration
S_G_corr_qRRHO = S_H_corr_qRRHO - T*S_S_tot_qRRHO

S_ept_qRRHO = e + S_E_tot_qRRHO
S_eph_qRRHO = e + S_H_corr_qRRHO
S_epg_qRRHO = e + S_G_corr_qRRHO

S_p_qRRHO = [N_epz,S_ept_qRRHO,S_eph_qRRHO,S_epg_qRRHO, S_S_tot_qRRHO]
S_p_qRRHO_corr = [S_H_corr_qRRHO, S_S_tot_qRRHO, S_G_corr_qRRHO]
#q_RRHO with 2D translational motion
#if args.mass1:
E_tot_qRRHO_2d = (E_v_qRRHO + E_t_2d)/4184 
H_corr_qRRHO_2d = E_tot_qRRHO_2d+ R*T/4184
S_tot_qRRHO_2d = (S_v_qRRHO + S_t_2d + R)/4184 #eneropy due to vibration and 2d translation
G_corr_qRRHO_2d = H_corr_qRRHO_2d - T*S_tot_qRRHO_2d
ept_qRRHO_2d = e + E_tot_qRRHO_2d
eph_qRRHO_2d = e + H_corr_qRRHO_2d
epg_qRRHO_2d = e + G_corr_qRRHO_2d
S_p_qRRHO_2d = [N_epz,ept_qRRHO_2d,eph_qRRHO_2d,epg_qRRHO_2d,S_tot_qRRHO_2d]
S_p_qRRHO_2dcorr = [H_corr_qRRHO_2d, S_tot_qRRHO_2d, G_corr_qRRHO_2d]
# print (epg_qRRHO_2d)

#q_RRHO with 1D translational motion
#if args.mass1:
E_tot_qRRHO_1d = (E_v_qRRHO + E_t_1d)/4184 
H_corr_qRRHO_1d = E_tot_qRRHO_1d+ R*T/4184
S_tot_qRRHO_1d = (S_v_qRRHO + S_t_1d + R)/4184 #eneropy due to vibration and 1d translation
G_corr_qRRHO_1d = H_corr_qRRHO_1d - T*S_tot_qRRHO_1d
ept_qRRHO_1d = e + E_tot_qRRHO_1d
eph_qRRHO_1d = e + H_corr_qRRHO_1d
epg_qRRHO_1d = e + G_corr_qRRHO_1d
S_p_qRRHO_1d = [N_epz,ept_qRRHO_1d,eph_qRRHO_1d,epg_qRRHO_1d,S_tot_qRRHO_1d]
S_p_qRRHO_1dcorr = [H_corr_qRRHO_1d, S_tot_qRRHO_1d, G_corr_qRRHO_1d]
 #save all the correction
S_q_RRHO_corr = [S_E_tot_qRRHO,S_H_corr_qRRHO,S_G_corr_qRRHO]   #Thermal, enthalpy and free energy correction for harmonic approximation with all degrees of freedom


#print S_p_qRRHO_corr
#print S_p_qRRHO_1dcorr
#print S_p_qRRHO_2dcorr
#
all_energy = [N_epz, S_epg, epg_1d, epg_2d, S_epg_qRRHO, epg_qRRHO_1d, epg_qRRHO_2d]
# electron energy plus zpe, free energy without qRRHO (immobile, 1D, 2D), free energy with qRRHO (immobile, 1D, 2D)

f.close()

# Write the new potential energy to txt file
name=str.split(file,'.')[0]+'_'+str(int(T))+'.txt'
f=open(name,'w')

f.write('The new temperature and pressure are %s K and %s Pa \n' % (T,P)) 
f.write('Total adsorbate: %s. The masses of the first and second adsorbate are %s amu and %s amu \n' % (n_ads, args.mass1, args.mass2)) 

 

f.write('# electron energy plus zpe, free energy without qRRHO (immobile, 1D, 2D), free energy with qRRHO (immobile, 1D, 2D):\n')
for item in all_energy:
    f.write(str(item)+'\n')

if args.mass1:
    f.write('If we consider vibrational degrees of freedom and 2D translational motion and use q_RRHO,the electronic and zero-point; electronic and thermal; electronic and enthalpy; electronic and free energy; S are:\n')
    for item in S_p_qRRHO_2d:
        f.write(str(item)+'\n')

if args.mass1:
    f.write('If we consider vibrational degrees of freedom and 1D translational motion and use q_RRHO,the electronic and zero-point; electronic and thermal; electronic and enthalpy; electronic and free energy; S are:\n')
    for item in S_p_qRRHO_1d:
        f.write(str(item)+'\n')


f.write('Use harmonic oscillator and keep all the degrees of freedom. The thermal energies (electronic, zpe, electronic and zero-point, electronic and thermal, electronic and enthalpy, electronic and free energy) in kcal/mol and S are:\n')
for item in N_p:
    f.write(str(item)+'\n')

f.write('Use harmonic oscillator and keep all the degrees of freedom. The thermal, enthalpy and free energy corrections are:\n')
for item in Harm_corr:
    f.write(str(item)+'\n')


#f.write('Use q_RRHO and keep all the degrees of freedom. The thermal energies (electronic and zero-point, electronic and thermal, electronic and enthalpy, electronic and free energy)and S in kcal/mol are:\n')
#for item in N_p_qRRHO:
#    f.write(str(item)+'\n')

f.write('If we only consider vibrational degrees of freedom and use harmonic oscillator, the electronic and zero-point, electronic and thermal, electronic and enthalpy, electronic and free energy are:\n')
for item in S_p:
    f.write(str(item)+'\n')

f.write('If we only consider vibrational degrees of freedom and use q_RRHO,the electronic and zero-point, electronic and thermal, electronic and enthalpy, electronic, free energy and S are:\n')
for item in S_p_qRRHO:
    f.write(str(item)+'\n')   

f.write('If we only consider vibrational degrees of freedom and use q_RRHO, The thermal, enthalpy and free energy corrections are:\n')
for item in q_RRHO_corr:
    f.write(str(item)+'\n')


   

f.close()
