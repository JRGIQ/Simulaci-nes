# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 08:52:29 2019

@author: Jheison Rene Gutierrez Gomez
"""
import matplotlib.pyplot as pt
import numpy as np
import math as mt

class isotermas():
    
    def puntos(self):
        self.P=np.array([105,204,300,401])
        self.Ve=np.array([1.96384,6.842,13.653,21.55])
        pt.plot(self.P,self.Ve,"k-")
        pt.grid(True)
        pt.xlabel("Presion")
        pt.ylabel("Volumen estandarizado")
        pt.show()
    def langmouir(self):
        self.XL=((1)/(self.P))
        self.YL=((1)/(self.Ve))
#        print("----------------------------------------------------------------")
#        print("las x,y para langmouir son : ",self.XL,self.YL)
#        
    def freundlich(self):
        self.XF=np.log(self.P)
        self.YF=np.log(self.Ve)
#        print("----------------------------------------------------------------")
#        print("las x,y para freundlich son : ",self.XF,self.YF)

        
    def BET(self):
        self.T=100
        self.A=13.8622 ### n-heptano
        self.B=2910.26
        self.C=216.432
        self.LnPsat=(self.A)-((self.B)/(self.T+self.C))
        self.Psat=np.exp(self.LnPsat)*((760)/(101.3*10))  #### para convertir la presion a cmHg
        self.XB=(self.P/self.Psat)
        self.YB=((self.P)/(self.Ve*(self.Psat-self.P)))
#        print("----------------------------------------------------------------")
#        print("las x,y para Bet son : ",self.XB,self.YB)
        
    def vector(self):
        self.X=self.XL
        self.Y=self.YL
        print("---------------------------------------------------------------")
        print("Los puntos ingresados son ejex-ejey respect/ : ",self.X,self.Y)
        
    def parametros_regresion(self):
        self.datos=int(input("ingrese el numero de datos : "))
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
        pt.xlabel("X")
        pt.ylabel("Y")
        pt.title("Regresion lineal")
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
        
        
isoterma=isotermas()
isoterma.puntos()
isoterma.langmouir()
isoterma.freundlich()
isoterma.BET()
isoterma.vector()
isoterma.parametros_regresion()
isoterma.coeficiente_regresion_R2()


