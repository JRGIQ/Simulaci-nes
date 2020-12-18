# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 21:49:01 2019

@author: JHEISON GUTIERREZ
"""


import numpy as np
import math as mt


class fugacidad():
    def parametros(self):
        self.R=8.314e-5 #m3*bar/K*mol
        self.T=473#float(input("ingrese la temperatura del sistema K : "))
        self.P=40#float(input("ingrese la presion del sistema bar : "))
        
        self.Tc1=304.2#float(input("ingrese la temperatura critica de la sustancia 1 (K) : "))
        self.Pc1=73.83#float(input("ingrese la presion critica de la sustancia 1 (bar) : "))
        self.w1=0.224#float(input("ingrese el factor acentrico de la sustancia 1"))
        self.Tr1=self.T/self.Tc1
        self.Pr1=self.P/self.Pc1
        
        self.Tc2=33.19#float(input("ingrese la temperatura critica de la sustancia 1 (K) : "))
        self.Pc2=13.13#float(input("ingrese la presion critica de la sustancia 1 (bar) : "))
        self.w2=-0.216#float(input("ingrese el factor acentrico de la sustancia 1"))
        self.Tr2=self.T/self.Tc2
        self.Pr2=self.P/self.Pc2
        
        self.Tc3=512.6#float(input("ingrese la temperatura critica de la sustancia 1 (K) : "))
        self.Pc3=80.97#float(input("ingrese la presion critica de la sustancia 1 (bar) : "))
        self.w3=32.042#float(input("ingrese el factor acentrico de la sustancia 1"))
        self.Tr3=self.T/self.Tc3
        self.Pr3=self.P/self.Pc3
        
        self.Tc4=647.1#float(input("ingrese la temperatura critica de la sustancia 1 (K) : "))
        self.Pc4=220.55#float(input("ingrese la presion critica de la sustancia 1 (bar) : "))
        self.w4=0.345#float(input("ingrese el factor acentrico de la sustancia 1"))
        self.Tr4=self.T/self.Tc4
        self.Pr4=self.P/self.Pc4
        
        self.Tc5=132.9#float(input("ingrese la temperatura critica de la sustancia 1 (K) : "))
        self.Pc5=34.99#float(input("ingrese la presion critica de la sustancia 1 (bar) : "))
        self.w5=0.048#float(input("ingrese el factor acentrico de la sustancia 1"))
        self.Tr5=self.T/self.Tc5
        self.Pr5=self.P/self.Pc5
        
        
    def coeficientes_viriales(self):
        self.B01=0.083-((0.422)/(self.Tr1**1.6))
        self.B11=0.139-((0.172)/(self.Tr1**4.2))
        self.B1=(self.R*self.Tc1/self.Pc1)*(self.B01+self.w1*self.B11)
        
        self.B02=0.083-((0.422)/(self.Tr2**1.6))
        self.B12=0.139-((0.172)/(self.Tr2**4.2))
        self.B2=(self.R*self.Tc2/self.Pc2)*(self.B02+self.w2*self.B12)
        
        self.B03=0.083-((0.422)/(self.Tr3**1.6))
        self.B13=0.139-((0.172)/(self.Tr3**4.2))
        self.B3=(self.R*self.Tc3/self.Pc3)*(self.B03+self.w3*self.B13)
        
        self.B04=0.083-((0.422)/(self.Tr4**1.6))
        self.B14=0.139-((0.172)/(self.Tr4**4.2))
        self.B4=(self.R*self.Tc4/self.Pc4)*(self.B04+self.w4*self.B14)
        
        self.B05=0.083-((0.422)/(self.Tr5**1.6))
        self.B15=0.139-((0.172)/(self.Tr5**4.2))
        self.B5=(self.R*self.Tc5/self.Pc5)*(self.B05+self.w5*self.B15)
        
    def coeficiente_fugacidad(self):
        self.LnΦ1=((self.B1*self.P)/(self.R*self.T))
        self.LnΦ2=((self.B2*self.P)/(self.R*self.T))
        self.LnΦ3=((self.B3*self.P)/(self.R*self.T))
        self.LnΦ4=((self.B4*self.P)/(self.R*self.T))
        self.LnΦ5=((self.B5*self.P)/(self.R*self.T))
        
        print("El coeficiente de fugacidad del CO2 es Φ1=",np.exp(self.LnΦ1))
        print("----------------------------------------------------------------")
        
        print("El coeficiente de fugacidad del H2 es Φ2=",np.exp(self.LnΦ2))
        print("----------------------------------------------------------------")
        
        print("El coeficiente de fugacidad del CH3OH es Φ3=",np.exp(self.LnΦ3))
        print("----------------------------------------------------------------")
        
        print("El coeficiente de fugacidad del H2O es Φ4=",np.exp(self.LnΦ4))
        print("----------------------------------------------------------------")
        
        print("El coeficiente de fugacidad del CO es Φ5=",np.exp(self.LnΦ5))
        print("----------------------------------------------------------------")
        
coeficiente=fugacidad()
coeficiente.parametros()
coeficiente.coeficientes_viriales()
coeficiente.coeficiente_fugacidad()
        
        
        
        
        
        
        
        
        
        
        
        