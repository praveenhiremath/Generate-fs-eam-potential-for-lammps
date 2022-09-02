# -*- coding: utf-8 -*-
from __future__ import division
import sys
import os
import subprocess
from subprocess import call
import numpy as np
import scipy 
import matplotlib
import pylab
from pylab import *
import math
from math import *
from decimal import Decimal
from eam_formalism import *

if os.path.exists("pair_interation_energies.dat"):
 os.remove("pair_interation_energies.dat")
if os.path.exists("local_density.dat"):
 os.remove("local_density.dat")
if os.path.exists("embed_energies.dat"):
 os.remove("embed_energies.dat")
if os.path.exists("APRIL13_Pure_W_fs.eam"):
 os.remove("APRIL13_Pure_W_fs.eam")
if os.path.exists("./dr_values.dat"):
 os.remove("./dr_values.dat")
if os.path.exists("drho_values.dat"):
 os.remove("drho_values.dat")


with np.errstate(divide='ignore'):
    np.float64(1.0) / 0.0

Nr=10000  ##Number of points at which electron density is evaluated
cutoff=6.000000000000   ##cutoff distance for all functions (Angstroms)

rc=3.25000000000000  # radial cutoff parameter

Nrho=Nr  ##Number of points at which the interatomic potential and embedding function is evaluated
dr=cutoff/Nr  #3.00000000000000E-0002   ##distance between points where the interatomic potential and embedding function is evaluated

print ("Nr = ",Nr," Nrho = ",Nrho," cutoff = ",cutoff)
print ("dr = ",dr)

'''
EAM parameters required to generate this potential are taken from the article --> G. Ackland, R. Thetford, An improved N-body semi-empirical model for
body-centred cubic transition metals, Phil. Mag. A 56 (1) (1987) 15–30.
'''
###parameters' values
A=1.896373
d=4.400224
c=3.25
c0=47.1346499
c1=-33.7665655
c2=6.2541999
B=90.3
b0=2.7411
alpha=1.2
num_species=1
at_nr=74
lat_par=3.1652   


###########APRIL 13 PHI CALCULATION
PHI=np.zeros(Nr)
pair_ene1=0.0
pair_ene2=0.0

for i in range(0,Nr,1):
  r=i*dr

eam_form=EAM(c,c0,c1,c2,B,b0,alpha,r)

for i in range(0,Nr,1):#  r=i*dr
  if ((r<b0) and (r<=rc)):  
    PHI[i]=eam_form.pair1(c,c0,c1,c2,B,b0,alpha,r)#*i*dr
  elif (r<=rc):
    PHI[i]=eam_form.pair2(c,c0,c1,c2,B,b0,alpha,r)
  else:
    PHI[i]=0.000000000000000E+00 

  
np.savetxt('pair_interation_energies.dat',PHI)


for i in range(0,Nr,1):
  PHI[i]=PHI[i]*i*dr

density=0.0
RHO=np.zeros(Nr)
for i in range(0,Nr,1):
  r=i*dr
  if (r<=d):
    RHO[i]=eam_form.local_density(d,r)
  else:
    RHO[i]=0.000000000000000E+00 

  
np.savetxt('local_density.dat',RHO)

rhomax=50.0
drho=rhomax/Nrho
print ("drho = ",drho)
F=np.zeros(Nrho)
emb_ene=0.0
for i in range(0,Nr,1):
  rho=i*drho
  F[i]=eam_form.embedd_energy(rho,A)
 
np.savetxt('embed_energies.dat',F)

##HEADER OF EAM/FS POTENTIAL FILE    
starting_lines=[None]*6
starting_lines[0]='#Potential parameters taken from Möller and Bitzek\n'
starting_lines[1]='#Interatomic potential for BCC Tungsten\n'
starting_lines[2]='#This file created by Praveenkumar Hiremath\n'
starting_lines[3]=str(num_species)+' W\n'       # '74 W\n'
starting_lines[4]=str(Nrho)+' '+str(drho)+' '+str(Nr)+' '+str(dr)+'0000000000 '+str(cutoff)+'\n'    #'Nrho, drho, Nr, dr, cutoff\n'
starting_lines[5]=str(at_nr)+' 183.84000000000000341061 3.16520000000000002842 BCC \n'

f1=open("APRIL13_Pure_W_fs.eam","w");
f1.writelines(starting_lines)
f1.close();

final_out='APRIL13_Pure_W_fs.eam'


f_handle = open(final_out, 'a')
np.savetxt(f_handle, F,fmt='%1.15E ')
f_handle.close()

f2_handle = open(final_out, 'a')
np.savetxt(f2_handle, RHO,fmt='%1.15E ')
f2_handle.close()

f3_handle = open(final_out, 'a')
np.savetxt(f3_handle, PHI,fmt='%1.15E ')
f3_handle.close()

r_ij=np.zeros(Nr)
for i in range(0,Nr,1):
  r_ij[i]=i*dr

np.savetxt('dr_values.dat',r_ij)

d_rho=np.zeros(Nrho)
for i in range(0,Nrho,1):
  d_rho[i]=i*drho

np.savetxt('drho_values.dat',d_rho)

print ("Potential file creation over! ")














