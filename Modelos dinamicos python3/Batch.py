# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 18:17:46 2020

@author: JRGIQ
"""

import numpy as np
import matplotlib.pyplot as pt

class bioreactor():
    
    
    def simulacion(self):
        
        
        
        self.K=0.025
        self.V=3.5
        self.muMax=0.178
        self.X0=7
        self.Y0=5
        self.t0=0
        self.t=15
        self.N=0.01
        self.h=int((self.t-self.t0)/self.N)
        self.tiempo=np.linspace(self.t0,self.t,self.h)
        self.concA=[]
        self.concB=[]


        for i in range(len(self.tiempo)):
            
    
            self.k=((self.muMax*self.X0)/(self.K+self.X0))
        
            self.dAdt=self.k*self.X0*self.V

            self.A=self.X0+self.dAdt*self.N        
                
            self.X0=self.A            
        
            self.concA.append(self.A*0.001)    


        pt.figure("BIOREACTOR BATCH")       
        pt.grid(True)   
        pt.plot(self.tiempo,self.concA,"k",label="Levadura",linewidth=2)
        pt.xlabel("Tiempo [h]")
        pt.ylabel("Masa [kg]")
        leg = pt.legend(loc='best', ncol=1, mode="center", shadow=True, fancybox=True)
        leg.get_frame().set_alpha(0.5)
    



simulacion=bioreactor()
simulacion.simulacion()



