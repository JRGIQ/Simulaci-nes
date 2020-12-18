# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 00:54:24 2019

@author: Jheison Gutierrez
"""

#### ESTE PROGRAMA CALCULA, GRAFICA: COORDENADA DE RxN VS ENERGIA DE GIBBS,CONSTANTE DE EQUILIBRIO


### PARA LA PRIMER ECUACION:     CO2 + 3H2 <---> CH3OH + H2O

##########################    A1=CO2   A2=H2   A3=CH3OH    A4=H2O

import math as mt
import numpy as np
from sympy import *
import matplotlib.pyplot as pt

#### Definir las variables simbolicas
x = Symbol('x')
y = Symbol('y')

#### se define el sistema de ecuaciones

#f=(4.0879*x**2)+(y**2)+(4.0879*x*y)+(0.2532*x)-(0.0439*y)-(0.1486)
#g=(1.12805*y**2)+(2.12805*x*y)-(0.12805*x)+(0.3048*y)-(0.4328)
f=(-16.29*x**4)-(3*x**3)+(4*x)+(22.06*x**3*y)-(y**4)-(5.6098*x**2*y)-(1.9219*x**2*y)+(2.5626*x*y**2)-(4*x*y)
g=(3*x**2)-(103259*y**2)-(206520*x)+(103265.84*x*y)-(103206)
####### se define la funcion matricial y la matriz
#### tambien la inversa del jacobiano
M=Matrix([f,g])
JI=(M.jacobian([x,y]))**-1


 ##### aproximacion inicial 
s=Matrix([0,1])
while(M.subs([(x,s[0]),(y,s[1])]).norm()>0.0001):
    
    s= (s )- (JI.subs([(x,s[0]),(y,s[1])])*M.subs([(x,s[0]),(y,s[1])]))
print (" x="+ str(s[0]))
print (" y = "+ str(s[1]))
print (M.subs([(x,s[0]),(y,s[1])]))#### matriz evaluada, tiene que ser cercana a cero.































