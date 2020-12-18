# -*- coding: utf-8 -*-
"""
Created on Sat Feb 16 19:16:05 2019

@author: Jheison Gutierrez
"""

from math import*
import numpy as np
import matplotlib.pyplot as pt

humedad=np.array([0.0,2.4,3.76,4.76,6.1,7.83,9.9,12.63,15.4,19.42]) ## para los dos lotes unidos
Pparcial=np.array([0.0,9.66,19.2,28.4,37.2,46.4,55,63.2,71.9,79.5])  ## para los 2 lotes unidos

Xw1=((humedad)/(100-humedad))
Yw1=((Pparcial)/(760-Pparcial))*((18.02)/(29))

pt.plot(Xw1,Yw1,"bo-")
pt.grid(True)
print("-------------------------------------------")
#print(Xw2,Yw2) ### para ver la tabla de puntos
pt.plot(0.2,0.0099,"ko")
pt.plot(0.149,0.0563,"ko")
a=[0.2,0.149]
b=[0.0099,0.0563]
pt.plot(a,b,"k")
coeficientes=np.polyfit(Yw1,Xw1,2)
polinomio=np.poly1d(coeficientes)
print("el polinomio Xw=f(Yw) es :",polinomio)