#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  2 15:31:36 2018

@author: jheison
"""

class calculo_volumen:
    
#    def_init_(self)
#        self.R=0
#        self.P=0
#        self.T=0
   
        
        
    def variables(self):
        self.R=83.15# cm3 bar mol-1 K-1
        self.T=float(input("ingrese temperatura (kelvin): "))
        self.P=float(input("ingrese presion (bar): "))
        print("________________________________________")
    
    def Volumen_ideal(self):
        self.Videal= self.R*self.T/self.P  
        print("el volumen ideal (cm3) es: ",self.Videal)
        print("_____________________________________________")
                                
        
    def Volumen_real(self):
        self.Vreal=((self.R*self.T)/(self.P))+self.B
        print("el valor del volumen real (cm3) es: ",self.Vreal)
        print("__________________________________________________")
    
    def volumen_real_ecuacion_cubica(self):
        self.V0=self.Videal  
        for iteration in range(1,100):
            self.V=((self.R*self.T)/(self.P))*(1+((self.B)/(self.V0))+((self.C)/(self.V0**2))) 
            self.V0=self.V
        print("el volumen real con la ecuacion cubica (cm3) es: ",self.V0)
        print("_____________________________________________________________")
        
    def segundo_coef_virial(self):
        self.B=-388
    
    def tercer_coef_virial(self):
        self.C=-26000
        
    def factor_compresibilidad(self):
        self.Z= self.Vreal/self.Videal
        print("el factor de compresibilidad es: ",self.Z)
        print("___________________________________________")
        
    def factor_compresibilidad_ec_estado(self):
        self.Zec=self.V0/self.Videal
        print("el factor de compresibilidad con ecuacion de estados es: ",self.Zec)
        print("_________________________________________________")
        
        
        
        
    
    
    
volumen1=calculo_volumen()
volumen1.variables()
volumen1.segundo_coef_virial()
volumen1.tercer_coef_virial()
volumen1.Volumen_ideal()
volumen1.Volumen_real()
volumen1.volumen_real_ecuacion_cubica()
volumen1.factor_compresibilidad()
volumen1.factor_compresibilidad_ec_estado()
