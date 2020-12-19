#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 14:40:05 2020

@author: jheison
"""

import numpy as np
import matplotlib.pyplot as pt

class bioreactor():
    
    def reactor(self):
                
        
        self.volumen=50#float(input("Volumen del reactor [m^3] = "))
        self.Y=0.43
        self.muMax= 0.178
        self.Ks=0.0025
        
        self.t0= 0
        
        self.t=500
        
        self.Sin=50
        self.F=3
        
        self.dt=0.001
        self.h=int((self.t-self.t0)/(self.dt))
        self.tiempo=np.linspace(self.t0,self.t,self.h)
        
        self.Xi=np.zeros(len(self.tiempo))
        self.Si=np.zeros(len(self.tiempo))
        self.mu=np.zeros(len(self.tiempo))

        
        
        self.S0=10
        self.X0=0.1
        
        self.X=self.X0
        self.S=self.S0
    
        
        
        for i in range(len(self.tiempo)):
            
            self.Xi[i]=self.X
            self.Si[i]=self.S
           
            
        #Ecuaciones constitutivas
        
            self.mu=self.muMax*((self.S0)/((self.S0)+self.Ks))

            
        # # Ecuaciones diferenciales
        
            self.dXdt=(self.mu*self.X0)-((self.X0*self.F)/(self.volumen))
            self.X=self.X0+self.dXdt*self.dt
            
            self.dSdt=((self.Sin*self.F)/(self.volumen))-((self.S0*self.F)/(self.volumen))-((self.mu*self.X0)/(self.Y))
            self.S=self.S0+self.dSdt*self.dt
            
            self.X0=self.X
            self.S0=self.S
            
    def graficas(self):
        
        pt.figure("Bioreactor CSTR",[5,4])
        pt.grid(True)
        pt.plot(self.tiempo,self.Xi,label="Levadura (Biomasa)",linewidth=2)
        pt.plot(self.tiempo,self.Si,label="Azúcares (Sustrato)",linewidth=2)
        pt.xlabel("Tiempo [h]")
        pt.ylabel("Concentración [kg/m^3]")
        pt.title("Volumen del reactor {} m^3".format(self.volumen))
        leg = pt.legend(loc='best', ncol=1, mode="center", shadow=True, fancybox=True)
        leg.get_frame().set_alpha(0.5)
        
                
        
        
simulacion=bioreactor()
simulacion.reactor()
simulacion.graficas()
        
 
    

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        