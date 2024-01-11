# -*- coding: utf-8 -*-
from __future__ import division
from math import *
'''
The code is for PhD research purposes

Author: Praveenkumar Hiremath
Email: praveenkumar.hiremath@mek.lth.se (Email at the University)
       praveenkumar.hiremath2911@gmail.com (Private email)
'''

class EAM:
    def __init__(self,c,c0,c1,c2,B,b0,alpha,r):
        self.__c_local = c
        self.__c0_local = c0
        self.__c1_local = c1
        self.__c2_local = c2
        self.__B_local = B
        self.__b0_local = b0
        self.__alpha_local = alpha 
        self.__r_local = r        #self.__pairene_local = pair_ene 

    
    def set(self,c,c0,c1,c2,B,b0,alpha,r):
        self.__c_local, self.__c0_local, self.__c1_local, self.__c2_local, self.__B_local, self.__b0_local, self.__alpha_local, self.__r_local = c,c0,c1,c2,B,b0,alpha,r  #, self.__pairene_local = pair_ene



###Function for pair potential
    def pair(self,c,c0,c1,c2,B,b0,alpha,r):
       constant=0.2407773468888542
       pair_ene=(((pow((r-c),2)*(c0+(c1*r)+c2*pow(r,2)))+(B*(pow((b0-r),3))*exp(-alpha*r)))+constant)#*neigh_atoms
       return pair_ene

    def pair1(self,c,c0,c1,c2,B,b0,alpha,r):
       constant1=0.0
       pair_ene1=(((pow((r-c),2)*(c0+(c1*r)+c2*pow(r,2)))+(B*(pow((b0-r),3))*exp(-alpha*r)))+constant1)
       return pair_ene1


    def pair2(self,c,c0,c1,c2,B,b0,alpha,r):
       constant2=0.0 
       pair_ene2=(((pow((r-c),2)*(c0+(c1*r)+c2*pow(r,2)))+(0.0*(pow((b0-r),3))*exp(-alpha*r)))+constant2)
       return pair_ene2


    def local_density(self,d,r):
       density=(pow((r-d),2))
       return density

    def embedd_energy(self,rho,A):
       emb_ene=-A*sqrt(rho)
       return emb_ene

