# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 22:42:56 2019

@author: Jheison Rene Gutierrez Gomez, JRIQ.
"""
import numpy as np
import matplotlib.pyplot as pt

### Tiempo de simulacion

class tanque_agitado():
    
    def parametros_entrada(self):
        
        self.a=0   ## Tiempo inicial minutos
        self.b=10  ## tiempo final minutos
        self.N=1000 ## Numero de iteraciones
        self.t=np.linspace(self.a,self.b,self.N) ## Crea un vector para el tiempo
        self.q=5 ## Flujo volumetrico de salida L/min
        self.qf=np.ones(len(self.t))*5.2 ## Flujo volumetrico del alimento L/min
        self.qf[400:]=5.1 ## Cambia el flujo del alimento en el segundo 4
        self.Caf=np.ones(len(self.t)) ## Concentracion del alimento
        self.Caf[500:]=0.5 ## Cambia la concentracion del alimento en el segundo 5
        self.Tf=np.ones(len(self.t))*300  ## Temperatura del alimento
        self.Tf[700:]=325 ## Cambia la temperatura del alimento en el segundo 7        
        
    def vectores(self):
        
        self.vectorCa=[] ## Crea un vector para las soluciones de la ecuacion diferecial dCa/dt
        self.vectorT=[]  ## Crea un vector para las soluciones de la ecuacion diferencial dT/dt
        self.vectorV=[]  ## Crea un vector para las soluciones de la ecuacion diferencial dV/dt
        
    def simulacion(self):
         
        self.V0=1 ## Volumen inicial de la mezcla
        self.Ca0=0  ## Concentracion inicial de la mezcla
        self.T0=500  ## Temperatura inicial de la mezcla
        self.rA=0  ## Constante de velocidad de RxN
        self.h=((self.b-self.a)/(self.N)) ## valor de una particion
       
        for i in range(len(self.t)):          
    
            self.dVdt=self.qf[i]-self.q ## Solucion ecuacion difereicial para el volumen
            self.dCadt=((self.qf[i]*self.Caf[i]-self.q*self.Ca0)/(self.V0))-self.rA-((self.Ca0*self.dVdt)/(self.V0)) ## Solucion ecuacion diferencial para la concentracion
            self.dTdt=((self.qf[i]*self.Tf[i]-self.q*self.T0)/(self.V0))-((self.T0*self.dVdt)/(self.V0)) ## Solucion ecuacion diferencial para la temperatura
            
            self.Cai=self.Ca0+self.dCadt*self.h ## Metodo euler
            self.Ti=self.T0+self.dTdt*self.h ## "  "
            self.Vi=self.V0+self.dVdt*self.h ## "  "
            self.Ca0=self.Cai ## Se renombra la concentracion en cada iteracion, solucion por pasos
            self.T0=self.Ti    ## " "
            self.V0=self.Vi  ## " "
            
            self.vectorCa.append(self.Ca0) ## Agrega cada solucion a un vector par la concentraicon
            self.vectorT.append(self.T0) ## Agrega cada solucion a un vector par la Temperatura
            self.vectorV.append(self.V0)   ## Agrega cada solucion a un vector par El volumen               
            
    def print(self): ## Modulo para graficar
        
        pt.figure(1,figsize=(9,9))
        
        pt.subplot(3,2,1)
        pt.plot(self.t,self.vectorT,"k-",linewidth=2)
        pt.xlabel("Tiempo min")
        pt.ylabel("Temperatura K")
        pt.grid(True)
        pt.legend(["Temperatura K"],loc='best')
        
        pt.subplot(3,2,2)
        pt.plot(self.t,self.Tf,"b-",linewidth=2)
        pt.xlabel("Tiempo min")
        pt.ylabel("Temperatura K, alimento")
        pt.grid(True)
        pt.legend(["Temperatura Tf"],loc='best')
        
        pt.subplot(3,2,3)
        pt.plot(self.t,self.vectorCa,"k-",linewidth=2)
        pt.xlabel("Tiempo min")
        pt.ylabel("Concentracion mol/L")
        pt.grid(True)
        pt.legend(["Concentracion Ca"],loc='best')
        
        pt.subplot(3,2,4)
        pt.plot(self.t,self.Caf,"b-",linewidth=2)
        pt.xlabel("Tiempo min")
        pt.ylabel("Concentracion mol/L, alimento")
        pt.grid(True)
        pt.legend(["Concentracion Caf"],loc='best')
        
        pt.subplot(3,2,5)
        pt.plot(self.t,self.vectorV,"k-",linewidth=2)
        pt.xlabel("Tiempo min")
        pt.ylabel("Volumen litros")
        pt.grid(True)
        pt.legend(["Volumen V"],loc='best')
        
        pt.subplot(3,2,6)
        pt.plot(self.t,self.qf,"b-",linewidth=2)
        pt.xlabel("Tiempo min")
        pt.ylabel("Velocidad de flujo L/min")
        pt.grid(True)
        pt.legend(["Velocidad de flujo Vf"],loc='best')
        
        
        pt.show()



tanque=tanque_agitado()
tanque.vectores()
tanque.parametros_entrada()
tanque.simulacion()
tanque.print()


