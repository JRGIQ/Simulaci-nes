# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 12:49:46 2019

@author: Jheison Rene Gutierrez Gomez, JRIQ.
"""

#### ESTE PROGRAMA CALCULA LA CONSTANTE DE EQUILIBRIO K, PARA UNA RxN HOMOGENEA BIMOLECULAR ejem 13.4 smith
import math as mt
import numpy as np
import matplotlib.pyplot as pt

class constante():
    def variables_del_sistema(self):
        self.T0=298.15#float(input("ingrese la temperatura inicial del sistema (K): "))
        self.T=500#np.array([298,398,498,598,698,798,898,989,1098,1198,1298,1398,1498]) #float(input("ingrese la temperatura final del sistema (K): "))
        self.tao=self.T/self.T0
        self.R=8.314
        
        
        self.A1=1.424#float(input("ingrese A del reactivo 1: "))
        self.B1=14.394*1e-03#float(input("ingrese B del reactivo 1 e-03 : "))*1e-03
        self.C1=-4.392*1e-06#float(input("ingrese C del reactivo 1 e-06 : "))*1e-06
        self.D1=0*1e+05#float(input("ingrese D del reactivo 1 e-05 : "))*1e+05
        self.ΔHo1=52510#float(input("ingrese el calor estandar de formacion del reactivo 1: "))
        self.ΔGo1=68460#float(input("ingrese la energia libre de Gibbs estandar (J/mol) del reactivo 1: "))
        
        self.A2=3.470#float(input("ingrese A del reactivo 2: "))
        self.B2=1.450*1e-03#float(input("ingrese B del reactivo 2 e-03 : "))*1e-03
        self.C2=0*1e-06#float(input("ingrese C del reactivo 2 e-06 : "))*1e-06
        self.D2=0.121*1e+05#float(input("ingrese D del reactivo 2 e-05 : "))*1e+05
        self.ΔHo2=-241818#float(input("ingrese el calor estandar de formacion del reactivo 2: "))
        self.ΔGo2=-228572#float(input("ingrese la energia libre de Gibbs estandar (J/mol) del reactivo 2: "))
        
        self.A3=3.518#float(input("ingrese A del producto 1: "))
        self.B3=20.001*1e-03#float(input("ingrese B del producto 1 e-03 : "))*1e-03
        self.C3=-6.002*1e-06#float(input("ingrese C del producto 1 e-06 : "))*1e-06
        self.D3=0*1e+05#float(input("ingrese D del producto 1 e-05 : "))*1e+05
        self.ΔHo3=-235100#float(input("ingrese el calor estandar de formacion del producto 1: "))
        self.ΔGo3=-168490#float(input("ingrese la energia libre de Gibbs estandar (J/mol) del producto 1: "))
        
        self.A4=0#float(input("ingrese A del producto 2: "))
        self.B4=0#float(input("ingrese B del producto 2 e-03 : "))*1e-03
        self.C4=0#float(input("ingrese C del producto 2 e-06 : "))*1e-06
        self.D4=0#float(input("ingrese D del producto 2 e-05 : "))*1e+05
        self.ΔHo4=0#float(input("ingrese el calor estandar de formacion del producto 2: "))
        self.ΔGo4=0#float(input("ingrese la energia libre de Gibbs estandar (J/mol) del producto 2: "))
        
        
        self.ΔA=(self.A3+self.A4)-(self.A1+self.A2)
        self.ΔB=(self.B3+self.B4)-(self.B1+self.B2)
        self.ΔC=(self.C3+self.C4)-(self.C1+self.C2)
        self.ΔD=(self.D3+self.D4)-(self.D1+self.D2)
        self.ΔHo298=(self.ΔHo3+self.ΔHo4)-(self.ΔHo1+self.ΔHo2)
        self.ΔGo298=(self.ΔGo3+self.ΔGo4)-(self.ΔGo1+self.ΔGo2)
                
        
    def integral_capacidad_calorifica_media(self): #### IDCPH smith van ness
        self.integralΔCp=((self.ΔA*self.T0)*(self.tao-1))+(self.ΔB*0.5*self.T0**2)*(self.tao**2-1)+(self.ΔC/3)*(self.T0**3)*(self.tao**3-1)+(self.ΔD/self.T0)*((self.tao-1)/(self.tao))
        self.integralΔCp2=self.ΔA*np.log(self.tao)+((self.ΔB*self.T0)+((self.ΔC*self.T0**2)+((self.ΔD)/(self.tao**2*self.T0**2)))*(self.tao+1)*(0.5))*(self.tao-1)
        print("el calor sensible medio ΔCp/R (J)  IDCPH es: ΔCp/R=",self.integralΔCp)
        print("----------------------------------------------------------------")
        print("el calor sensible medio ΔCp/RT  IDCPS (J) es: ΔCp°/RT=",self.integralΔCp2)
        print("----------------------------------------------------------------")
        print("La entalpia estandar de reaccion a 298K (J/mol) es: ΔHo298=",self.ΔHo298)
        print("----------------------------------------------------------------")
        print("La energia libre de Gibbs estandar de reaccion a 298K es (J/mol) es: ΔG°298=",self.ΔGo298)
        print("----------------------------------------------------------------")
        
    def energia_gibbs_estandar_RxN(self):
        self.ΔGo=((self.ΔGo298-self.ΔHo298)/(self.R*self.T0))+((self.ΔHo298)/(self.R*self.T))+((1/self.T)*(self.integralΔCp))-(self.integralΔCp2)
        print("La energia libre de Gibbs estandar de reaccion es: ΔG°/RT=",self.ΔGo)
        print("----------------------------------------------------------------")
        
    def constante_de_equilibrio(self):
        self.LnK=-(self.ΔGo)
        self.LnKo298=-((self.ΔGo298)/(self.R*298.15))
        self.K=np.exp(self.LnK)
        self.Ko298=np.exp(self.LnKo298)
        print("La constante de equilibrio a temperatura T es: K=",self.K)
        print("----------------------------------------------------------------")
        print("La constante de equilibrio a condiciones estandar es: K°298=",self.Ko298)
        print("----------------------------------------------------------------")
        
    def grafica_ΔGo_vs_T(self):
        
        pt.xlabel('Temperatura (K)')
        pt.ylabel('ΔG° (J/mol)')
        pt.grid(True)
        pt.title(' Temperatura vs ΔG° RxN')
        pt.plot(self.T,self.ΔGo*self.R*self.T,'k-',linewidth=2)
        pt.show()
        
    def grafica_LnK_vs_T(self):
        pt.xlabel('Temperatura (K) 1/T')
        pt.ylabel('constante de equilibrio LnK')
        pt.grid(True)
        pt.title('Temperatura vs Ln(K)')
        pt.plot(1/self.T,self.LnK,'b-',linewidth=2)
        pt.show()
        
    def grafica_K_vs_T(self):
        pt.xlabel('Temperatura (K) T')
        pt.ylabel('constante de equilibrio K')
        pt.grid(True)
        pt.title('Temperatura vs constante de equilibrio K')
        pt.plot(self.T,self.K,'r-',linewidth=2)
        pt.show()
    def grafica_H_vs_T(self):
        pt.xlabel('Temperatura (K) T')
        pt.ylabel('Entalpia de reaccion')
        pt.grid(True)
        pt.title('Temperatura vs Entalpia')
        pt.plot(self.integralΔCp,self.T,'b-',linewidth=2)
        pt.show()
        



constante=constante()
constante.variables_del_sistema()
constante.integral_capacidad_calorifica_media()
constante.energia_gibbs_estandar_RxN()
constante.constante_de_equilibrio()
#constante.grafica_ΔGo_vs_T()
#constante.grafica_LnK_vs_T()
#constante.grafica_K_vs_T()
#constante.grafica_H_vs_T()
        
        
        