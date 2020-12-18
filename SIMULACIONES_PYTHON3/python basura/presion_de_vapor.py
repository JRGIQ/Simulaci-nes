#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 01:51:25 2018

@author: jheison
"""

#def presion_de_vapor():
#    valores={}
#    continuar= "s"
#    while continuar=="s":
#        temperatura=int(input("ingrese la temperatura"))
#        presion=int(input("ingrese la presion"))
#        valores[temperatura]=presion
#        continuar=input("quiere agregar otro valor?: s/n")
#    return valores
#valores=presion_de_vapor()
#
##print(diccionario())
#
#def consultar(valores):
#    temp=input("ingrese la temperatura")
#    if temp in valores:
#        print("la presion de vapor es: ",valores[temp])

from control.matlab import*
import matplotlib.pyplot as plt

g=tf(1,[1,1,1])
print(g)
t,y=step_response(g)
plt.plot(t,y)
plt.grid()
plt.xlabel('tiempo')
plt.ylabel('respuesta escalon')