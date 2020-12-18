#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 09:56:37 2018

@author: jheison
"""



#import numpy as np
#def funcion(x):
#    y= x**2
#    return y
#x=np.array([1,2,3,4,5,6,7,8,9,10])
#print(funcion(x))


#Calculo del volumen de un cilindro V(r,h)
#import math
#def funcion(r,h):
#    V=math.pi*r*r*h
#    return V
#r= int(input("ingrese el valor del radio: "))
#h=int(input("ingrese el valor de la altura: "))
#volumen_cilindro=funcion(r,h)
#print("el volumen del cilindro es: ",volumen_cilindro)



##########################################################################
#####DEFINIR FUNCION DEFINIENDO VARIABLES COMO LISTAS, Y RECORRERLAS
#def sumar(x,y,*lista):
#    suma=x+y
#    for i in range(len(lista)):
#        suma=suma+lista[i]
#    return suma
#print(sumar(1,2,3))

##########:############################################################
###DEFINIR UN ELEMENTO EN ESPECIAL DE LA LISTA(MAYOR)
#def mayor(lista):
#    may=lista[0]
#    for i in range(5,len(lista)): 
#        if lista[i]>may:
#            may=lista[i]
#    return may
#lista=[10,2,3,4,56,6,7,8,9]
#print("el numero mayor es: ", mayor(lista))
###################################################################
###AÑADIR ELEMENTOS A TUPLAS....EN UNA LISTA
#def diccionario():
#    diccionary={}
#    continuar= "s"
#    while continuar=="s":
#        español=input("ingrese la palabra en español")
#        ingles=input("ingrese la palabra en ingles")
#        diccionary[español]=ingles
#        continuar=input("quiere agregar otra palabra s/n")
#    return diccionary
#
#print(diccionario())

#def consultar(diccionary):
#    palabra=input("ingrese palabra en español")
#    if palabra in diccionary:
#        print("la palabra en ingles es: ",diccionary[palabra])



#
#from control import*
#from control.matlab import*
#import math as mt
#from scipy import* 
#import matplotlib.pyplot as plt
#g=tf(1,[1,1,1])
#print(g)
#t,y=step_response(g)
#plt.plot(t,y)
#plt.grid()
#plt.xlabel('tiempo')
#plt.ylabel('respuesta escalon')
#
#        
        
      ##################para cuadros de texto en la grafica  
        
#pt.xlabel('presion')
#pt.ylabel('factor de compresibilidad')
#pt.grid(True)
#leg = pt.legend(loc='lower right', ncol=1, mode="center", shadow=True, fancybox=True)
#leg.get_frame().set_alpha(0.7)
#pt.show()
#coeficientes=np.polyfit(presion,factorZ_360K,3)
#polinomio=np.poly1d(coeficientes)
#print(polinomio)
#x=Symbol('x')
#print(integrate(-3.744e-05*x**3 + 0.0004824*x**2 - 0.01677*x + 0.9979,x))
#f=lambda x:-9.36e-6*x**4 + 0.0001608*x**3 - 0.008385*x**2 + 0.9979*x
#print(quad(f,0,15.41))
#        best
#        upper right
#        upper left
#        lower left
#        lower right
#        right
#        center left
#        center right
#        lower center
#        upper center
#        center

#---------------------------------------------------------------------------
#coeficientes=np.polyfit(,factorZ_340K,3)
#polinomio=np.poly1d(coeficientes)
#print(polinomio)