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







################SOLUCION DEL CURSO CON ODEINT
#import numpy as np
#import matplotlib.pyplot as plt
#from scipy.integrate import odeint
#
## define mixing model
#def vessel(x,t,q,qf,Caf,Tf):
#    # Inputs (4):
#    # qf  = Inlet Volumetric Flowrate (L/min)
#    # q   = Outlet Volumetric Flowrate (L/min)
#    # Caf = Feed Concentration (mol/L)
#    # Tf  = Feed Temperature (K)
#
#    # States (3):
#    # Volume (L)
#    V = x[0]
#    # Concentration of A (mol/L)
#    Ca = x[1]
#    # Temperature (K)
#    T = x[2]
#
#    # Parameters:
#    # Reaction
#    rA = 0.0
#
#    # Mass balance: volume derivative
#    dVdt = qf - q
#
#    # Species balance: concentration derivative
#    # Chain rule: d(V*Ca)/dt = Ca * dV/dt + V * dCa/dt
#    dCadt = (qf*Caf - q*Ca)/V - rA - (Ca*dVdt/V)
#
#    # Energy balance: temperature derivative
#    # Chain rule: d(V*T)/dt = T * dV/dt + V * dT/dt
#    dTdt = (qf*Tf - q*T)/V - (T*dVdt/V)
#
#    # Return derivatives
#    return [dVdt,dCadt,dTdt]
#
## Initial Conditions for the States
#V0 = 1.0
#Ca0 = 0.0
#T0 = 350.0
#y0 = [V0,Ca0,T0]
#
## Time Interval (min)
#t = np.linspace(0,10,100)
#
## Inlet Volumetric Flowrate (L/min)
#qf = np.ones(len(t))* 5.2
#qf[50:] = 5.1
#
## Outlet Volumetric Flowrate (L/min)
#q = np.ones(len(t))*5.0
#
## Feed Concentration (mol/L)
#Caf = np.ones(len(t))*1.0
#Caf[30:] = 0.5
#
## Feed Temperature (K)
#Tf = np.ones(len(t))*300.0
#Tf[70:] = 325.0
#
## Storage for results
#V  = np.ones(len(t))*V0
#Ca = np.ones(len(t))*Ca0
#T  = np.ones(len(t))*T0
#
## Loop through each time step
#for i in range(len(t)-1):
#    # Simulate
#    inputs = (q[i],qf[i],Caf[i],Tf[i])
#    ts = [t[i],t[i+1]]
#    y = odeint(vessel,y0,ts,args=inputs)
#    # Store results
#    V[i+1]  = y[-1][0]
#    Ca[i+1] = y[-1][1]
#    T[i+1]  = y[-1][2]
#    # Adjust initial condition for next loop
#    y0 = y[-1]
#
## Construct results and save data file
#data = np.vstack((t,qf,q,Tf,Caf,V,Ca,T)) # vertical stack
#data = data.T             # transpose data
#np.savetxt('data.txt',data,delimiter=',')
#
## Plot the inputs and results
#plt.figure()
#
#plt.subplot(3,2,1)
#plt.plot(t,qf,'b--',linewidth=3)
#plt.plot(t,q,'b:',linewidth=3)
#plt.ylabel('Flow Rates (L/min)')
#plt.legend(['Inlet','Outlet'],loc='best')
#
#plt.subplot(3,2,3)
#plt.plot(t,Caf,'r--',linewidth=3)
#plt.ylabel('Caf (mol/L)')
#plt.legend(['Feed Concentration'],loc='best')
#
#plt.subplot(3,2,5)
#plt.plot(t,Tf,'k--',linewidth=3)
#plt.ylabel('Tf (K)')
#plt.legend(['Feed Temperature'],loc='best')
#plt.xlabel('Time (min)')
#
#plt.subplot(3,2,2)
#plt.plot(t,V,'b-',linewidth=3)
#plt.ylabel('Volume (L)')
#plt.legend(['Volume'],loc='best')
#
#plt.subplot(3,2,4)
#plt.plot(t,Ca,'r-',linewidth=3)
#plt.ylabel('Ca (mol/L)')
#plt.legend(['Concentration'],loc='best')
#
#plt.subplot(3,2,6)
#plt.plot(t,T,'k-',linewidth=3)
#plt.ylabel('T (K)')
#plt.legend(['Temperature'],loc='best')
#plt.xlabel('Time (min)')
#
#plt.show()








