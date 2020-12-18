# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 17:59:15 2019

@author: Jheison Rene Gutierrez Gomez, "JRIQ".
"""
import numpy as np
import matplotlib.pyplot as pt


class regresion_lineal():
    
    def puntos(self):
        self.datos=int(input("Ingrese el numero de datos n : "))
        self.vectorx=[]
        self.vectory=[]

        for i in range(self.datos):
            self.x=float(input("Ingrese el dato para eje x : "))
            self.y=float(input("Ingrese el dato para eje y : "))
            self.vectorx.append(self.x)
            self.vectory.append(self.y)
#            
    def vector(self):
        self.X=np.array(self.vectorx)
        self.Y=np.array(self.vectory)
        print("---------------------------------------------------------------")
        print("Los puntos ingresados son ejex-ejey respect/ : ",self.X,self.Y)
        
    def parametros_regresion(self):
        self.sumax=sum(self.X)
        self.sumay=sum(self.Y)
        self.xcuadrado=sum(self.X**2)
        self.ycuadrado=sum(self.Y**2)
        self.sumaxy=sum(self.X*self.Y)
        self.promx=sum(((self.X)/(self.datos)))
        self.promy=sum(((self.Y)/(self.datos)))
        self.pendiente=(self.sumax*self.sumay-self.datos*self.sumaxy)/(self.sumax**2-self.datos*self.xcuadrado)
        self.intercepto=self.promy-self.pendiente*self.promx
        print("---------------------------------------------------------------")
        print("La pendiente de la recta que se ajusta a los datos es : ",self.pendiente)
        print("---------------------------------------------------------------")
        print("El intercepto de la recta que se ajusta a los datos es : ",self.intercepto)
        pt.plot(self.X,self.pendiente*self.X+self.intercepto,label="Ajuste con regresion lineal")
        pt.xlabel("tiempo (min)")
        pt.ylabel("mg/L Aminas aromaticas")
        pt.title("Cinetica segundo orden")
        pt.plot(self.X,self.Y,"o",label="Datos experimentales")
        pt.plot(label="R^2 = 0.9999")
        pt.grid(True)
        leg=pt.legend(loc="upper left",ncol=1,mode="center", shadow=True)
        leg.get_frame().set_alpha(0.5)
        
    def coeficiente_regresion_R2(self):
        self.sigmax=np.sqrt(self.xcuadrado/self.datos-self.promx**2)
        self.sigmay=np.sqrt(self.ycuadrado/self.datos-self.promy**2)
        self.sigmaxy=self.sumaxy/self.datos-self.promx*self.promy
        self.R2=(self.sigmaxy/(self.sigmax*self.sigmay))**2
        print("---------------------------------------------------------------")
        print("El coeficiente de regresion lineal R^2 es : ",self.R2)
        print("---------------------------------------------------------------")
          
            

regresion=regresion_lineal()
regresion.puntos()
regresion.vector()
regresion.parametros_regresion()
regresion.coeficiente_regresion_R2()

    
