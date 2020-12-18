# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 18:50:35 2019

@author: Jheison Rene Gutierrez Gomez."JRIQ".
"""
import matplotlib.pyplot as pt

   ##### Ejemplo TdeMasa ejercicio 2.3 notas de clase cap2
#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 1-  PARAMETROS Y CONSTANTES.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 1.1-  PARAMETROS DEL EQUIPO.

Ma=14 ##Masa del soluto A
areaM=140 ## Area de la TdeMasa
ME=20   ##Masa inicial Kg de (Extracto+soluto)
MR=5  ##Masa inicial Kg de (Refinado+soluto)
Y0=0.05  ## concentracion inicial del soluto en el refinado Gas
X0=0.0095  ## Concentracion inicial del soluto en el Extracto Liquido
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 1.2-  PARAMETROS DE LAS SUSTANCIAS Y CONSTANTES.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------


WRs=((1)/(1-Y0)) # Concentracion de refinado (GAS) libre de soluto en la corriente de entrada
WEs=((1)/(1+X0)) # Concentracion de extracto (LIQUIDO) libre de soluto en la corriente de entrada
Kxw=1.3e-5#float(input("Ingrese el coeficiente de TdeMasa local englobante Kxw: "))
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 4-  CREACIÓN DE VECTORES PARA GUARDAR DATOS.


vector1=[]  ## Crea un vector vacio para guardar los datos iterados para el tiempo
vector2=[]  ## Crea un vector vacio para guardar los datos iterados para el soluto en solvente puro (Kg /Kg Liquido puro) 
vector3=[]  ## Crea un vector vacio para guardar los datos iterados para la cantidad de (extracto+soluto)
vector4=[]  ## Crea un vector vacio para guardar los datos iterados para la cantidad de (refinado+soluto)
vector5=[]  ## Crea un vector vacio para guardar los datos iterados para el soluto en solvente puro (Kg /Kg Gas puro)
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 3-  TIEMPO DE SIMULACIÓN, PASO Y NÚMERO DE ITERACIONES.


t1=0 ## tiempo de inicio
t2=int(input("Ingrese el tiempo de duracion de la sumulacion (segundos) : ")) ## Pide el tiempo de simulacion final
N=t2 ##numero de particiones N=t2 para que el tiempo sea de a un segundo
h=((t2-t1)/(N)) #### numero de intervalos para el metodo euler



w1=MR  ## Kg iniciales del Refinado para el metodo euler
w12=ME ##Kg iniciales del Extracto para el metodo euler



for i in range(t1,t2): ### comienza el ciclo for, todo lo que ingrese a este ciclo se itera. 
    Mtrans=Kxw*(0.036-1.25*X0)*Ma*areaM  ## calcula la TdeMasa, tambien se recalcula en cada iteracion 
    Xi=((Mtrans+X0*w1)/(Mtrans+w1)) ## Recalcula la concentracion del soluto en el solvente liquido
    Yi=Y0-((Mtrans)/(w1)) ##Recalcula la concentracion del soluto en el refinado
    wi=w1-(h*Mtrans)  ### Aplicacion del metodo euler para el refinado
    wi2=w12+(h*Mtrans) ##Aplicacion del metodo euler para el extracto
    
    w1=wi  ## la cantidad de (refinado+soluto) que se calcula hasta este paso, se convierte en el dato de entrada de la siguiente iteracion.
    X0=Xi  ## La concentracion de soluto en el extraxto Liquido puro que se calcula hasta este paso,se convierte en el dato de entrada de la siguiente iteracion.
    Y0=Yi  ## La concentracion de soluto en el refinado Gas puro que se calcula hasta este paso,se convierte en el dato de entrada de la siguiente iteracion.
    w12=wi2  ##la cantidad de (Extracto+soluto) que se calcula hasta este paso, se convierte en el dato de entrada de la siguiente iteracion.
         
        
    tiempo=t1+i*h    ### Calcula el tiempo en cada iteracion
    vector1.append(tiempo) ### Guarda los valores de cada iteracion en un vector para poderlas graficar
    vector2.append(X0)   ### "     "
    vector3.append(w1)  ### "     "
    vector4.append(w12)  ### "     "
    vector5.append(Y0)  ### "     "
    
    ### Grafica para la TdeMasa concentracion del soluto Xw-Yw en el solvente puro, a travez del tiempo.
pt.figure(1,figsize=(7,5)) ## llama una figura, y le pone medidas.
pt.plot(vector1,vector2,'r-',linewidth=2,label="Kg soluto A / Kg solvente liquido puro   Xw")## Grafica el vector concentracion en la fase extracto (liquida), pone color y forma a la linea
pt.grid(True)  ## grafica la cuadrilla 
pt.xlabel('Segundos') ##pone titulo en la coordenada del eje x
pt.ylabel('Xw , Yw')  ##pone titulo en la coordenada del eje y
pt.title('Simulacion concentracion de soluto A en el solvente puro')  ## pone titulo a la grafica 
pt.subplot(111)## comando para dibujar en la misma grafica, poniendo tamaño a su vez.
pt.plot(vector1,vector5,'b-',linewidth=2,label="Kg soluto A / Kg solvente Gas puro Yw ") ## Grafica el vector concentracion en la fase refinado (gas), pone color y forma a la linea
pt.grid(True) ### pone cuadrilla a la grafica
leg = pt.legend(loc='lower right', ncol=1, mode="center", shadow=True, fancybox=True) ## pone condiciones para las leyendas en el grafico.
leg.get_frame().set_alpha(0.7)
pt.show() ### grafica una unica imagen hasta aca.

###  Grafica la TdeMasa soluto+refinado en funcion del tiempo
pt.figure(2,figsize=(7,5))#### llama una figura, y le pone medidas.
pt.grid(True)## pone cuadrilla a la grafica
pt.plot(vector1,vector3,'y-',linewidth=2,label="Kg de mezcla (soluto+refinado) Rws")## Grafica el vector  refinado+soluto (gas), pone color y forma a la linea
pt.xlabel('Segundos') ## pone nombre al eje x
pt.ylabel('Rws') ## pone nombre al eje y
pt.title('Simulacion TdeMasa Kg mezcla (Refinado+soluto) Rws') ## pone titulo a la grafica
leg = pt.legend(loc='upper right', ncol=1, mode="center", shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.7)

    ###  Grafica la TdeMasa soluto+extracto en funcion del tiempo  
pt.figure(3,figsize=(7,5))
pt.grid(True)
pt.plot(vector1,vector4,'g-',linewidth=2,label="Kg de mezcla (soluto+Extracto) Ews")## Grafica el vector concentracion en la extracto+soluto (liquido), pone color y forma a la linea
pt.xlabel('Segundos')
pt.ylabel('Ews')
pt.title('Simulacion TdeMasa Kg mezcla (Extracto+soluto) Ews')
leg = pt.legend(loc='lower right', ncol=1, mode="center", shadow=True, fancybox=True)
leg.get_frame().set_alpha(0.7)
pt.show()




        
        
            
            
            

        
        
        
        
   
#        
        