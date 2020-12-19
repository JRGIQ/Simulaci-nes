# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 21:02:37 2020

@author: JRGIQ
"""

#  MODULO 0-  LIBRERÍAS DE PYTHON.

import numpy as np ### Libreria de python que permite realizar operaciones con vectores.
import matplotlib.pyplot as pt# Libreria de python que permite realizar gráficas.

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 1-  PARAMETROS Y CONSTANTES.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 1.1-  PARAMETROS Y DEL EQUIPO.

Volumen = 0.1*0.065*0.006  ##  Volumen del IdeCalor [m^3].
SP = 37+273.15  ##  Set Point [K].

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 1.2-  PARAMETROS DE LAS SUSTANCIAS Y CONSTANTES.

Densidad=1.192  ##  Densidad aire[kg/m^3].
Masa=Volumen*Densidad  ##  Masa de agua dentro de la coraza [kg].
Cp_L=910 ##   Cp de la placa [Kj/Kg K].

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 2-  CONDICIONES INICIALES. (Estado estacionario).(ini=Inicial, in=Entrada, out=Salida).

Pert_ini = 50  ##  Apertura incial de la válvula de entrada de líquido (d) [%].
T_in = 25+273.15  ##  Temperatura de entrada (agua) al tanque [k].
Val_Ctrl_ini = 50  ## Apertura inicial de la válvula de vapor [%].
T_out_ini=SP  ##  Temperatura inicial de salida del proceso.
FlujoM_aire_in=2.7472527e-7 ##  Flujo másico de líquido que ingresa al IdeCalor [kg/s].
Q_0=0.003 # Calor de la placa en estado estacionario [Kw]
#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 3-  TIEMPO DE SIMULACIÓN, PASO Y NÚMERO DE ITERACIONES.

t_Inicial = 0  ##  Tiempo en que inicia la simulación [s].
t_Final =2000#int(input("Ingrese el tiempo de simulación [s] = ")) ## Tiempo en que finaliza la simulación [s].
Paso = 1 ##float(input("ingrese el tamaño del paso para el método Euler [s] = "))  ## Tamaño de cada partición de intervalo de tiempo.
N_Part=int((t_Final-t_Inicial)/(Paso))  ##  Determina el número de particiones necesarias para resolver el método.
tiempo = np.linspace(t_Inicial, t_Final, N_Part)  ##  Crea vector de tiempo para poder graficarlo [s].
Iter = N_Part  ##  Número de iteraciones, se puede decir que es igual al número de particiones.
t_Pert = int(0.5 * N_Part)  ##  Tiempo y valor en que se aplica la perturbación [%].
n=1#int(input("Ingrese el tiempo de muestreo del controlador"))

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 4-  CREACIÓN DE VECTORES PARA GUARDAR DATOS.

T_out_i = np.zeros(len(tiempo))  ##  Crea un vector de datos para la temperatura de salida, variable de proceso (y).
Pert_i = np.zeros(len(tiempo))  ##  Crea un vector de datos para la (acción) o entrada de aire (d).
Error_i = np.zeros(len(tiempo))  ##  Crea un vector de datos para el error= set point- T_out.
Val_Ctrl_i = np.zeros(len(tiempo)+1)  ##  Crea un vector de datos para la (acción) de control o apertura del válvula de vapor (u).
SP_i=np.zeros(len(tiempo))  ##  Crea un vector de datos para el Set Point.
FlujoM_aire=np.zeros(len(tiempo))  ##  Crea un vector de datos para el flujo de aire (Fluido de proceso). 
Q_i=np.zeros(len(tiempo))  ##  Crea un vector de datos para el calor energía.

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 5-  SET POIN, PERTURBACIONES, Y VARIABLE DE CONTROL.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 5.1-  SET POINT.

SP_i[t_Inicial:]=SP  ##  Llena el vector del Set Point con el valor asignado como referencia.

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 5.2-  PERTURBACIONES.

Pert_i[t_Inicial:] = Pert_ini  ##  Llena el vector para la perturbación con el valor inicial asignado.
Pert_i[t_Pert:] = Pert_ini + 30  ##  Cambia el valor de la perturbación al vector a un tiempo t_Pert(tiempo perturbación)

#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 6-  PARÁMETROS DEL CONTROLADOR.

Kp=-100 ##  Parámetro de ganancia proporcional [%].
ti=5  ##  Parámetro para el tiempo integral [s].
td=2  ##  Parámetro para el tiempo derivativo [s].
ITC=1  ##  Intervalo de tiempo de la acción de control.
e_1k=0  ##  Error en el instante de tiempo anterior (i-1).
e_2k=0  ##  Error en el instante de tiempo anterior (i-2).

#________________________________________________________________________________________________________________________________________________________________________
## AQUÍ EMPIEZA LA SOLUCIÓN ITERATIVA DEL MODELO DEL PROCESO CON SU CONTROLADOR.##

T_out= T_out_ini  ##  Asigna el valor inicial a la temperatura para iniciar el método [K].

for i in range(0, Iter, 1):
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 7-  ACTUALIZACIÓN DE LAS VARIABLES DE PROCESO, CONTROL Y ERROR.

    T_out_i[i] = T_out  ##  Guarda los datos para cada solución de la Ecuación diferencial de Temperatura.
    Val_Ctrl_i[t_Inicial] = Val_Ctrl_ini  ##  Guarda los datos para cada solución de la variable de control (Apertura de la válvula de vapor).
    Error_i[i] = (SP_i[i] - T_out_i[i])  ##  Guarda el error en cada iteración.
  
#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 8-  SOLUCIÓN DE ECUACIONES.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 8.1-  SOLUCIÓN DE ECUACIONES CONSTITUTIVAS.
    
    FlujoM_aire[i] = (Val_Ctrl_i[i] / 100) *2* FlujoM_aire_in  ##  Soluciona y actualiza la ecuación para el flujo de proceso (perturbación).    
    Q_i[i]=(Pert_i[i]/100)*2*Q_0 ## Actualiza la perturbación.
 
    
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 8.2-  SOLUCIÓN DE ECUACIONES DIFERENCIALES APLICANDO EL MÉTODO EULER.
    
    dTdt =((FlujoM_aire[i]/Masa)*(T_in-T_out_i[i]))+((1)/(Masa*Cp_L))*Q_i[i]  ##  Soluciona la derivada de la ED para la temperatura.
    T_out= T_out_ini + dTdt * Paso  ##  Soluciona la ecuación diferencial, agregandole a la temperatura inicial el cambio que tiene esta en una pequeña partición o instante de tiempo (Paso).
    T_out_ini = T_out  ##  Actualiza la temperatura para iniciar con la siguiente partición. 
    
#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 9-  SELECCIÓN Y ACCIONES DEL CONTROLADOR.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 9.1-  CONTROLADOR (P).

    if (td==0 and ti>=9999):  ##  Control (P), depende de los parametros td y ti.
        Delta_u = Kp*ITC*Error_i[i]  ##  Ecuación para el control (P).

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 9.2-  CONTROLADOR (PI),(PD) Y (PID).

    else:                
         Delta_u= Kp*ITC*(Error_i[i]-e_1k) + Kp*(ITC/ti)*e_1k + Kp*(td/ITC)*(Error_i[i]-2*e_1k+e_2k)  ##  Ecuación para el control (PID)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 9.3-  ACTUALIZACION DEL ERROR.
         
    e_1k=Error_i[i]  ##  Actualiza el error anterior (i-1).
    e_2k=e_1k  ##  Actualiza el error anterior (i-2).
        
#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 9.4-  ACCIÓNES DE CONTROL.    
    
    if (Val_Ctrl_i[i]+Delta_u)>100: 
        Val_Ctrl_i[i+1]=100  ##  Corrige la acción de control para que no pase de su límite máximo.
        
    elif (Val_Ctrl_i[i]+Delta_u)<0:        
        Val_Ctrl_i[i+1]=0  ##  Corrige la acción de control para que no pase de su límite mínimo.
    
    else:
         Val_Ctrl_i[i+1]=Val_Ctrl_i[i]+Delta_u  ##  Acción de control dentro de los límites establecidos, se actualiza normalmente sumandole el cambio.  

Val_Ctrl_i=Val_Ctrl_i[:i+1]  ##  Ajusta el vector a las dimensiones requeridas. 
      
#________________________________________________________________________________________________________________________________________________________________________
#  MODULO 10-  GRAFICAS Y RESULTADOS.        
pt.figure("MODELO INTERCAMBIADOR DE CALOR (CONTROL PID)", [10,5])                 
           
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.2-  GRAFICA LOS DATOS, TEMPERATURA DE SALIDA (y) (Variable de proceso) VS TIEMPO .
                                                                                                          #_______#
                                                                                                          #|1 2 3|#
pt.subplot(2, 3, 1)  ##  Hace una sub-figura (N°filas , N°columnas , Posicion de la figura en el subplot) #|4 5 6|#.
pt.tight_layout(pad=4, w_pad=5, h_pad=2)  ##  Espaciado entre las figuras del subplot (Espaciado entre la margen , Espaciado entre figuras (Horizontal) , Espaciado entre figuras (Vertical).
pt.plot(tiempo, T_out_i-273.15, "k", linewidth=2)  ##  Grafica los datos requeridos (x , y , Color de la linea , Grosor de la linea).
pt.plot(tiempo, SP_i-273.15, "y--", label="Set point", linewidth=1)  ##  Grafica el Set Point.
pt.grid(True)  ##  Agrega la cuadrilla a la grafica.
pt.xlabel("Tiempo [s]")  ##  Agrega título al eje x.
pt.ylabel("Temperatura\nfinal [°C] ")  ##  Agrega título al eje y.
pt.ylim(min(T_out_i-273.15)-0.1,max(T_out_i-273.15)+0.1)  ##  Ajusta los limites en el eje (y) para mejor visualización.
pt.title("Variable de proceso")
pt.tick_params(labelsize=8)  ##  Ajusta el tamaño de los títulos para los ejes  (x , y), y de igual manera ajusta el tamaño de los números en cuadrilla(Grid).
pt.legend(loc="best",prop={'size': 7})  ##  Agrega la leyenda en la grafica.("best" es la mejor ubicación).

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.3-  GRAFICA LOS DATOS, ACCIÓN DE PERTURBACIÓN (d) VS TIEMPO .

pt.subplot(2, 3, 2) 
pt.plot(tiempo, Pert_i, "r", linewidth=2)
pt.grid(True) 
pt.xlabel("Tiempo [s]") 
pt.ylabel("Calor [%]")
pt.title("Perturbación") 
pt.ylim(min(Pert_i)-10,max(Pert_i)+10)  
pt.tick_params(labelsize=8)  

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.4-  GRAFICA LOS DATOS, ACCIÓN DE CONTROL (u) VS TIEMPO .

pt.subplot(2, 3, 3)
pt.plot(tiempo, Val_Ctrl_i, "b", linewidth=2) 
pt.grid(True)
pt.xlabel("Tiempo [s]")
pt.ylabel("Flujo aire \nventilador  [%]") 
pt.ylim(min(Val_Ctrl_i)-10,max(Val_Ctrl_i)+10) 
pt.title("Variable de control")
pt.tick_params(labelsize=8)  

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.5-  GRAFICA LOS DATOS, POTENCIA O CALOR TRANSFERIDO VS TIEMPO.

pt.subplot(2, 3, 4)
pt.plot(tiempo, Error_i, "c", linewidth=2) 
pt.grid(True) 
pt.xlabel("Tiempo [s]") 
pt.ylabel("Error [°C]") 
pt.ylim(min(Error_i)-0.01,max(Error_i)+0.01)  
pt.tick_params(labelsize=8)  

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.6-  GRAFICA LOS DATOS, FLUJO FLUIDO DE PROCESO (Variable de perturbación) VS TIEMPO.

pt.subplot(2, 3, 5)
pt.plot(tiempo, FlujoM_aire, "y", linewidth=2)
pt.grid(True)
pt.xlabel("Tiempo [s]")
pt.ylabel("Flujo aire [kg/s]") 
pt.tick_params(labelsize=8)  

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  MODULO 10.7-  GRAFICA LOS DATOS, FLUJO DE VAPOR (Variable Ctrl) VS TIEMPO.

pt.subplot(2, 3, 6)
pt.plot(tiempo, Q_i, "g", linewidth=2) 
pt.grid(True) 
pt.xlabel("Tiempo [s]") 
pt.ylabel("Flujo de calor [kw]") 
pt.ylim(min(Q_i)-0.01,max(Q_i)+0.01)  
pt.tick_params(labelsize=8) 

#__________________________________________________________________________________________________________________#________________________________________________________________________________________________________________________________________________________________________
print(Masa)