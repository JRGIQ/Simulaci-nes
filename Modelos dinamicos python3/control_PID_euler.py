# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 13:35:49 2019

@author: Jheison Rene Gutierrez Gomez, JRIQ.
"""

import numpy as np
import matplotlib.pyplot as pt
from scipy.integrate import odeint

class controlador_PI():
    def parametros(self):
        
        ### ecuacion diferencial "taop*dy/dt=-y+Kp(u(t)-thetap)"
        self.Kp=2  
        self.taop=200  
        self.thetap=0
        self.a=0 ## Tiempo inicial segundos
        self.b=1200  ## Tiempo final segundos
        self.N=1200  ## Numero de particiones 
        self.t=np.linspace(self.a,self.b,self.N) ## Crea vector para el tiempo
    
    def matrices(self):
        self.vp=np.zeros(len(self.t)) ## Crea vector para la variable de proceso
        self.vp1=np.zeros(len(self.t)) ## Crea vector para la solucion de la ecuacion diferencial
        self.Set=np.zeros(len(self.t)) ## crea vector para el set point
        self.Set[50:500]=10 ## valor del set point despues de 50 segundos
        self.Set[500:]=5## Cambia el valor del set point despues de 500 segundos
        self.error=np.zeros(len(self.t)) ## Crea vector para el error
        self.integral=np.zeros(len(self.t)) ## Crea vector para el termino integral del PI
        self.derivativo=np.zeros(len(self.t)) ## Crea vector para el termino derivativo para el PID
        self.cp=np.zeros(len(self.t)) ## Crea vector para la variable de control
        self.P=np.zeros(len(self.t)) ## Crea vector para el termino proporcional
        self.I=np.zeros(len(self.t)) ## Crea vector para el termino integral
        self.D=np.zeros(len(self.t)) ## Crea vector para el termino derivativo
        
        
    def PID(self):
        self.Kc=2 ## Constante proporcional
        self.taoI=10 ## Constante integral--con 100 no hay oscilaciones
        self.taoD=0 ## Constante derivativa
        
    def simulacion_modo_auto(self):
        self.y0=0 ## Valor inicial PVI con control
        self.y01=0 ## Valor inicial PVI sin control
        self.deltat=self.t[1]-self.t[0] ## Delta t, para el termino integral y derivativo
        self.h=((self.b-self.a)/(self.N)) ## numero de particiones del intervalo para el metodo euler
        self.cp_hi=100 ## Valor maximo de la variable de control 100 %
        self.cp_lo=0 ## Valor minimo de la variable de control 0 %
        
        for i in range(len(self.t)):
            
            self.dydt1=((-self.y01)/(self.taop))+(self.Kp/self.taop) ## Ecuacion diferencial 
            self.y1=self.y01+(self.dydt1*self.h) ## Metodo euler sin control
            self.vp1[i]=self.y1 ## Ingresa los resultados al vector creado
            self.y01=self.y1 ## Renombra la entrada para la nueva iteracion
            
            
            self.error[i]=self.Set[i-1]-self.vp[i-1]  ## Error entre el set point y la salida          
            self.integral[i] =( self.integral[i-1] + self.error[i])* self.deltat ## Termino integral
            self.derivativo[i]=((self.vp[i]-self.vp[i-1])/(self.deltat)) ## Termino derivativo
            
            self.P[i]= self.Kc*self.error[i] ## Guarda el termino proporcional en el vector creado
            self.I[i]=((self.Kc)/(self.taoI))*self.integral[i] ## Guarda el termino integral en el vector creado
            self.D[i]=-self.Kc*self.taoD*self.derivativo[i] ## Guarda el termino derivativo en el vector creado
            self.cp[i] = self.cp[0] + self.P[i] + self.I[i] + self.D[i] ## Guarda la variable de control en el vector creado
            
            self.dydt=((-self.y0)/(self.taop))+(self.Kp/self.taop)*self.cp[i] ## Ecuacion diferencial con control
            self.y=self.y0+(self.dydt*self.h) ## Metodo euler con control
            self.y0=self.y ## Renombra la entrada para la nueva iteracion
            self.vp[i]=self.y ## Guarda los datos de la variable de proceso en el vector creado       
                        
            
            if self.cp[i] > self.cp_hi: ## Limite superior
                self.cp[i] = self.cp_hi
                self.integral[i] = (self.integral[i] - self.error[i])* self.deltat # anti-reset windup
                
            if self.cp[i] < self.cp_lo:  ## limite inferior
                self.cp[i] = self.cp_lo
                self.integral[i] = (self.integral[i] - self.error[i] )* self.deltat # anti-reset windup
            ## Condiciona los limites de la variable de proceso
            

                
                
    def graficas(self): ## Modulo creado unicamente para graficar
        pt.figure(1,figsize=(5,10))
        
        pt.subplot(3,1,1) ### Para la variable de proceso
        pt.grid(True)
        pt.xlabel("Tiempo (segundos)")
        pt.ylabel("Variable de proceso")
        pt.plot(self.t,self.vp,"k-",label="VP")
        pt.plot(self.t,self.Set,"b--",label="Set point")
        pt.legend(loc="best")
        
        pt.subplot(3,1,2) ## Para la variable de control
        pt.grid(True)
        pt.xlabel("Tiempo (segundos)")
        pt.ylabel("Variable controlada")
        pt.plot(self.t,self.cp,"r-",label="VC")
        pt.legend(loc="best")
        
        pt.subplot(3,1,3) ## Para la respuesta a lazo abierto
        pt.grid(True)
        pt.xlabel("Tiempo (segundos)")
        pt.ylabel("Salida")
        pt.plot(self.t,self.vp1,"g-",label="Respuesta lazo abierto")
        pt.legend(loc="best")           
        print(len(self.cp))
        
control=controlador_PI()
control.parametros()
control.matrices()
control.PID()
control.simulacion_modo_auto()
control.graficas()    