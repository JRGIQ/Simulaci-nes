# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 16:42:04 2019

@author: Jheison Rene Gutierrez Gomez, JRIQ.
"""

##### CONTROLADOR P, PARA EL NIVEL DE UN TANQUE

import numpy as np
import matplotlib.pyplot as pt
from scipy.integrate import odeint

class control_P_tanque():
    def parametros(self):
        self.Ro=1000 ## Densidad Kg/m3
        self.A=1  ## Area transversal del tanque m
        self.C=50  ## Coeficiente de la valvula Cv
        self.Ubios=0
        self.Kc=50 ## Constante controlador proporcional
        self.a=0  ## Tiempo de inicio seg
        self.b=10  ## Tiempo final seg
        self.N=100  ## Numero de particiones del intervalo
        self.SP=10
        self.t=np.linspace(self.a,self.b,self.N) ## Vector para el tiempo
        
    def matrices(self):
        self.V=np.zeros(len(self.t)) ## Matriz para la valvula
        self.E=np.zeros(len(self.t)) ## Matriz para el error Sp-Valor variable salida
        self.L=np.zeros(len(self.t)) ## matriz para en nivel del fluido
        self.Set=np.zeros(len(self.t)) ## Matriz para el set point
        self.Set[60:]=self.SP ## Inicia el set point en 6 segundos
        
    def simulacion(self):
        
        self.h=((self.b-self.a)/(self.N)) ## numero de particiones para el metodo euler
        self.nivel0=0 ## Nivel de inicio del fluido
#        self.valvula=0
#        self.valvula0=0
#        self.SP=10 ## Set point del controlador metros

        for i in range(len(self.t)-1):
            
            self.error=self.Set[i]-self.nivel0 ## Error 
            self.valvula=self.Ubios+self.Kc*self.error ## Ecuacion controlador P, para la valvula
            self.valvula=max(0,self.valvula) ## Define la apertura minima de la valula %
            self.valvula=min(100,self.valvula) ## Define la apertura maxima de la valvula %
            self.dhdt=((self.C)/(self.Ro*self.A))*self.valvula ## Ecuacion diferencial para la altura
            self.nivel=self.nivel0+self.dhdt*self.h ## Metodo euler para el nivel del fluido
            self.nivel0=self.nivel ## Renombra la entrada para la nueva iteracion
            
            self.L[i+1]=self.nivel0 ## Graba el resultado de cada iteracion en la matriz para el nivel      
            self.E[i+1]=self.error  ## Graba el resultado de cada iteracion en la matriz para el error
            self.V[i+1]=self.valvula  ## Graba el resultado de cada iteracion en la matriz para la respuesta de la valvula
        
        pt.grid(True) ## Agrega cuadrilla a la grafica    
        pt.xlabel("Tiempo (segundos)") ## Titulo para el eje x
        pt.ylabel("Nivel (metros)") ## Titulo para el eje y
        pt.plot(self.t,self.L,"b--",label="Nivel") ## grafica tiempo vs nivel
        pt.plot(self.t,self.Set,"k-",label="Set point") ## Grafica tiempo vs set point
        pt.legend(loc="best") ## Pone el cuadro de informacion en la grafica.
         

tanque=control_P_tanque()
tanque.parametros()
tanque.matrices()
tanque.simulacion()