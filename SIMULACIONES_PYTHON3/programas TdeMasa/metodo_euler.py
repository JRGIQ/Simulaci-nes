# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 22:55:42 2019

@author: Jheison Gutierrez
"""      
    
    
#euler(0,1200,5,funcion,1200)

### CALCULO DEL POLINOMIO Yw=f(Xw) CON LA DCE PARA HALLAR LOS VALORES LIMITES

#x=np.array([0,0.25,0.5,1.25,2.5])*1e-3
#y=np.array([0,2.8,6.2,19,54])*1e-3
#coeficientes=np.polyfit(x,y,2)
#polinomio=np.poly1d(coeficientes)
#pt.grid(True)
#pt.plot(x,y,'ko-',linewidth = 1,label="Yw vs Xw")
#leg = pt.legend(loc='upper left', ncol=1, mode="center", shadow=True, fancybox=True)
#leg.get_frame().set_alpha(0.7) 
#print(polinomio)
#pt.show()

#EJEMPLO 5.2 TREYBAL

humedad=np.array([0.0,2.4,3.76,4.76,6.1,7.83,9.9,12.63,15.4,19.42]) ## para los dos lotes unidos
Pparcial=np.array([0.0,9.66,19.2,28.4,37.2,46.4,55,63.2,71.9,79.5])  ## para los 2 lotes unidos
#humedadjabon=np.array([0.0,2.4,3.76,4.76,6.1])
#Pparcial=np.array([0.0,9.66,19.2,28.4,37.2])
#humedadjabon2=np.array([7.83,9.9,12.63,15.4,19.42])
#Pparcial2=np.array([46.4,55,63.2,71.9,79.5])

#Yw1=np.array((Pparcial/760))
#Y2w=np.array((Pparcial2/760))

Xw1=((humedad)/(100-humedad))
#Xw2=((humedadjabon2)/(100-humedadjabon2))
Yw1=((Pparcial)/(760-Pparcial))*((18.02)/(29))
#Yw2=((Pparcial2)/(760-Pparcial2))*((18.02)/(29))
pt.plot(Xw1,Yw1,"bo-")
#pt.plot(Xw2,Yw2)
pt.grid(True)
print("-------------------------------------------")
#print(Xw2,Yw2) ### para ver la tabla de puntos
pt.plot(0.2,0.0099,"ko")
pt.plot(0.149,0.0563,"ko")
a=[0.2,0.149]
b=[0.0099,0.0563]
pt.plot(a,b,"k")
coeficientes=np.polyfit(Yw1,Xw1,2)
coeficientes2=np.polyfit(Xw1,Yw1,2)
#coeficientes2=np.polyfit(Yw2,Xw2,2)
polinomio=np.poly1d(coeficientes)
polinomio2=np.poly1d(coeficientes2)
print("el polinomio Xw=f(Yw) es :",polinomio)
print("el polinomio Yw=X(Yw) es :",polinomio2)












