# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 23:34:48 2019

@author: Jheison Rene Gutierrez Gomez,JRIQ.
"""



import matplotlib.pyplot as pt ## Libreria para realizar graficas
import numpy as np ## Libreria para realizar operaciones con matrices
import sys ## Libreria para que el codigo corra en algun programa online


class Simulador_etilbenceno():###orden sustancias=[etileno,etano,propileno,benceno,tolueno,etilbenceno,dietilbenceno]
    
    def parametros(self):
        
        print("***********************************************************************")
        print("[P E N G-- R O B I N S O N]")
        print("***********************************************************************")
        self.TC=80  ### Temperatura de proceso °C
        self.T=self.TC+273.15 ### Temperatura K
        self.P=210 ### Presion de proceso Kpa
        self.R=8.314  ### Constante de los gases J/mol K
        self.M=np.array([28.05,30.07,42.08,78.11,92.14,106.17,134.22]) ### g/mol          
        self.w=np.array([0.08722,0.09469,0.14332,0.20989,0.26497,0.30348,0.40257])
        self.Tc=np.array([282.3,305.3,365.6,562.2,591.8,617.2,657.87])##K
        self.Pc=np.array([50.40,48.72,46.65,48.98,41.06,36.06,28.3])*100##Kpa
        self.fracciones_mol=np.array([0.0018644293973655234, 0.024137240210021574, 0.006896354345720449, 0.6132452071199341, 1e-9, 0.3181102058113509, 0.03574656311560738])### Fracciones molares a la entrada [vapor]
        self.fracciones_mol_Liquido=np.array([0,0,0,0,0,0,0])    ### Fracciones molares a la entrada [liquido]
        self.corriente_alimento= 290 #### Kmol/h
        
        
        ##### CONSTANTES DE ANTOINE pag 235, propierties gases and liquids, 5ta edicion 
        
        self.KA=np.array([3.91382,3.95405,3.95606,3.98523,4.05043,4.06861,4.12598])
        self.KB=np.array([596.5260,663.720,789.6240,1184.240,1327.62,1415.77,1592.590])
        self.KC=np.array([256.370,256.681,247.580,217.572,217.62500,212.300,202.440])

        ####### Creacion de vectores para llenarlos posteriormente
        
        self.fw=np.ones(len(self.fracciones_mol))
        self.Tr=np.ones(len(self.fracciones_mol))
        self.Pr=np.ones(len(self.fracciones_mol))
        self.alpha=np.ones(len(self.fracciones_mol))
        self.a_pura=np.ones(len(self.fracciones_mol))
        self.b_pura=np.ones(len(self.fracciones_mol))
        self.A_pura=np.ones(len(self.fracciones_mol))
        self.am=np.ones((len(self.fracciones_mol),len(self.fracciones_mol)))
        self.amL=np.ones((len(self.fracciones_mol),len(self.fracciones_mol)))
        self.bmL=np.ones(len(self.fracciones_mol))
        self.bm=np.ones(len(self.fracciones_mol))
        self.wm=np.ones(len(self.fracciones_mol))
        
        
    def parametros_Ecuacion_Estado(self):
        
        self.Tr=self.T/self.Tc  ### Temperatura reducida
        self.fw=0.37464+1.5422*self.w-0.26992*self.w**2  ### Funcion fw Peng-Robinson Tabla 2.5 [Seader]
        self.alpha=(1+self.fw*(1-self.Tr**0.5))**2  
        #self.alpha[posicionH2]=1.202*np.exp(-0.30288*self.Tr[3])######## Utilizar unicamente para el hidrogeno en la mezcla
        
        self.a_pura=(0.457248*((self.R*self.Tc))**2/(self.Pc))*self.alpha ### Funcion a Peng-Robinson Tabla 2.5 [Seader]
        self.b_pura=0.07780*((self.R*self.Tc)/(self.Pc))  ### Funcion a Peng-Robinson Tabla 2.5 [Seader]
        
        for i in range(len(self.fracciones_mol)):
            
            self.bm[i]=self.fracciones_mol[i]*self.b_pura[i]  ## Regla de mezcla bi [vapor] Ecuacion 2-50 Seader 
            self.b_mezcla=sum(self.bm) ## Regla de mezcla bi [vapor] Ecuacion 2-50 Seader
            
            self.bmL[i]=self.fracciones_mol_Liquido[i]*self.b_pura[i] ## Regla de mezcla bi [liquido] Ecuacion 2-50 Seader 
            self.b_mezcla_L=sum(self.bmL) ## Regla de mezcla bi [liquido] Ecuacion 2-50 Seader 
            

            for j in range(len(self.fracciones_mol)):
                
                self.am[i][j]=self.fracciones_mol[i]*self.fracciones_mol[j]*(self.a_pura[i]*self.a_pura[j])**0.5 ## Regla de mezcla aij [vapor] Ecuacion 2-49 Seader 
                self.amL[i][j]=self.fracciones_mol_Liquido[i]*self.fracciones_mol_Liquido[j]*(self.a_pura[i]*self.a_pura[j])**0.5 ## Regla de mezcla aij [vapor] Ecuacion 2-49 Seader 
        self.a_mezcla=sum(sum(self.am)) 
        self.a_mezcla_L=sum(sum(self.amL)) ###Realiza la doble sumatoria de la ecuaicon 2-49 [Seader]
        
    def factor_compresibilidad(self):
        
        self.A_mezcla=((self.a_mezcla*self.P)/(self.R*self.T)**2)  ### Funcion A mezcla vapor, ecuacion 2-47 [Seader]
        self.A_mezcla_L=((self.a_mezcla_L*self.P)/(self.R*self.T)**2) ### Funcion A mezcla liquida, ecuacion 2-47 [Seader]
        
        self.B_mezcla=((self.b_mezcla*self.P)/(self.R*self.T))### Funcion B mezcla vapor, ecuacion 2-48 [Seader]
        self.B_mezcla_L=((self.b_mezcla_L*self.P)/(self.R*self.T))### Funcion B mezcla liquida, ecuacion 2-48 [Seader]
        
        self.Ai=((self.a_pura*self.P)/(self.R*self.T)**2)### Funcion A pura, ecuacion 2-47 [Seader]
        self.Bi=((self.b_pura*self.P)/(self.R*self.T))### Funcion B pura, ecuacion 2-48 [Seader]
        
        self.coeficiente_Z3=1  ### Coeficiente de Z3
        self.coeficiente_Z2=-(1-self.B_mezcla)### Coeficiente de Z2
        self.coeficiente_Z=self.A_mezcla-3*self.B_mezcla**2-2*self.B_mezcla### Coeficiente de Z
        self.coeficiente=-(self.A_mezcla*self.B_mezcla-self.B_mezcla**2-self.B_mezcla**3)### Coeficiente constante
        
        self.coeficientes=[self.coeficiente_Z3,self.coeficiente_Z2,self.coeficiente_Z,self.coeficiente] ### Vector que une los coeficientes Z
        self.raices=np.roots(self.coeficientes) #### Funcion que calcula las raices de la ecuacion cubica [Z1 Z2 Z3]

        
    def raices_reales_imaginarias(self):
        
        if self.raices[0].imag==self.raices[1].imag==self.raices[2].imag==0:  ### Condicion unica si todas las raices del polinomio son reales
            print("_________________________________________________________")
            print("######TODAS RAICES REALES######")
            print("_________________________________________________________")

            self.Z_vapor=max(self.raices.real)  ### Escoge la raiz mayor para el vapor
            self.Z_liquido=min(self.raices.real) ### Escoge la raiz menor para el liquido
            
            
        elif self.raices[0].imag==0: ### Condicion unica si solo la primer raiz posee el valor real
            
            print("_________________________________________________________")
            print("######UNA SOLA RAIZ REAL######")
            print("_________________________________________________________")
            self.Z_vapor=self.raices[0]  ### Escoge el primer valor como la raiz de vapor
            self.Z_liquido=0 ### Condicion si no hay raiz de liquido
            
            
        elif self.raices[1].imag==0: ### Condicion unica si solo la segunda raiz posee el valor real
            
            print("_________________________________________________________")
            print("######UNA SOLA RAIZ REAL######")
            print("_________________________________________________________")
            self.Z_vapor=self.raices[1] ### Escoge el primer valor como la raiz de vapor
            self.Z_liquido=0 ### Condicion si no hay raiz de liquido
            
        elif self.raices[2].imag==0:### Condicion unica si solo la tercer raiz posee el valor real
            
            print("_________________________________________________________")
            print("######UNA SOLA RAIZ REAL######")
            print("_________________________________________________________")
            self.Z_vapor=self.raices[2] ### Escoge el primer valor como la raiz de vapor
            self.Z_liquido=0 ### Condicion si no hay raiz de liquido
            
    def coeficiente_fugacidad(self):           
            
            
            for k in range(len(self.fracciones_mol_Liquido)):
                           
                if self.fracciones_mol_Liquido[k]==0:  ### Condicion para una unica raiz real [Si no ingresa liquido]
                    
                    self.Plog=self.KA-((self.KB)/(self.TC+self.KC)) ### Presion de saturacion [Antoine]
                    self.Psat=(10**self.Plog)*100 ## Presin de saturacion Kpa
                    
                    self.coeficientes_K=(self.Psat/self.P) ### Coeficientes K para condicion supercritica o cuando no ingresa vapor al equipo. Tabla 2, ecuacion 3 [Seader]
                    
                    
                       
                elif self.fracciones_mol_Liquido[k]!=0:  #### Condicion para calculo de valores K en funcion de la relacion de fugacidades---- Cuando todas las raices son reales y ingresa liquido al equipo                  
                    
                    self.Ln_fi=(self.Z_vapor-1)*(self.Bi/self.B_mezcla)-np.log(self.Z_vapor-self.B_mezcla)-(self.A_mezcla/self.B_mezcla)*(2*(self.Ai/self.A_mezcla)**0.5)*np.log(1+(self.B_mezcla/self.Z_vapor)) ### Calcula el coeficiente de fugacidad de cada componente en la fase vapor
                    self.fi=np.exp(self.Ln_fi) ### Calcula el coeficiente de fugacidad de cada componente en la fase vapor
            
                    self.Ln_fiL=(self.Z_liquido-1)*(self.Bi/self.B_mezcla_L)-np.log(self.Z_liquido-self.B_mezcla_L)-(self.A_mezcla_L/self.B_mezcla_L)*(2*(self.Ai/self.A_mezcla_L)**0.5)*np.log(1+(self.B_mezcla_L/self.Z_liquido)) ### Calcula el coeficiente de fugacidad de cada componente en la fase liquida
                    self.fiL=np.exp(self.Ln_fiL) ### Calcula el coeficiente de fugacidad de cada componente en la fase liquida
                    self.coeficientes_K=self.fiL/self.fi ### Calcula los coeficientes K en funcion de la relacion de coe_fugacidad_liquido/coe_fugacidad_vapor Ecuacion 2-26 [Seader]
                    
    def Rachford_Rice(self):   ### ALGORITMO SEPARADOR FLASH--- Rachford Rice---Tabla 4.4 Seader
        
        self.iteraciones=100  ### Define numero de iteraciones
        self.nv0=0  ### Valor inicial de la iteracion [0-1]
        
        for z in range(self.iteraciones):  #### METODO DE NEWTON RAPSON PARA RESOLVER EL METODO DE RACHFORD-RICE
        

            self.f_newton=((self.fracciones_mol*(1-self.coeficientes_K))/(1+(self.nv0*(self.coeficientes_K-1)))) ### Funcion f metodo Newton 
            self.f_suma=sum(self.f_newton) ### Sumatoria de f
    
            self.f_prima_newton=(((self.fracciones_mol)*((1-self.coeficientes_K)**2))/(((self.nv0*(self.coeficientes_K-1))+1)**2)) ### Funcion f'[derivada f] metodo Newton
            self.f_prima_suma=sum(self.f_prima_newton) ### Sumatoria de f
            
            self.nv_newton=self.nv0-((self.f_suma)/(self.f_prima_suma)) ### Metodo de Newton-Rapson
            self.nv0=self.nv_newton ## Se renombra el parametro para nueva iteracion hasta convergencia
                
            if self.nv0<=0: ### Condicion para que no aparezca error de division por cero
                self.nv0=1e-9
            elif self.nv0>1:  ### Condicion cuando nvo converge a un valor mayor que uno... recordando que el valor maximo es uno.
                self.nv0=1
            
                
            self.tolerancia=abs((self.nv_newton-self.nv0)/(self.nv0))  ### Tolerancia metodo Newton
                
        if self.tolerancia<=0.001:  ### Se le da un valor a la tolerancia deseada
            
            print("___________________________________________________________________")
            print("La fraccion de vapor ********[RACHFORD-RICE]******* [nv] es : ",self.nv0) 
            print("___________________________________________________________________")
            
            
    def fracciones_salida_Flash(self):
        
        self.Flujo_vapor=self.corriente_alimento*self.nv0  ## Define el flujo de vapor en la corriente de salida
        self.Flujo_liquido=self.corriente_alimento-self.Flujo_vapor  ## Define el flujo de liquido en la corriente de salida
        
        self.fraccion_mol_vapor=((self.fracciones_mol*self.coeficientes_K)/(1+self.nv0*(self.coeficientes_K-1))) ### Define las fracciones de cada componente en la fase vapor
        self.fraccion_mol_liquido=((self.fracciones_mol)/(1+self.nv0*(self.coeficientes_K-1))) ### Define las fracciones de cada componente en la fase liquida
        
        
        #### SE DEFINEN LAS FRACCIONES MOLARES PARA CADA UNO DE LOS COMPONENTES DE LA MEZCLA A LA SALIDA DEL SEPARADOR FLASH
        
        self.fraccion_mol_etileno_Fase_liquida=round(self.fraccion_mol_liquido[0],4)
        self.fraccion_mol_etileno_Fase_vapor=round(self.fraccion_mol_vapor[0],4)
        
        self.fraccion_mol_etano_Fase_liquida=round(self.fraccion_mol_liquido[1],4)
        self.fraccion_mol_etano_Fase_vapor=round(self.fraccion_mol_vapor[1],4)
        
        self.fraccion_mol_propileno_Fase_liquida=round(self.fraccion_mol_liquido[2],4)
        self.fraccion_mol_propileno_Fase_vapor=round(self.fraccion_mol_vapor[2],4)
        
        self.fraccion_mol_benceno_Fase_liquida=round(self.fraccion_mol_liquido[3],4)
        self.fraccion_mol_benceno_Fase_vapor=round(self.fraccion_mol_vapor[3],4)
        
        self.fraccion_mol_tolueno_Fase_liquida=round(self.fraccion_mol_liquido[4],4)
        self.fraccion_mol_tolueno_Fase_vapor=round(self.fraccion_mol_vapor[4],4)
        
        self.fraccion_mol_etilbenceno_Fase_liquida=round(self.fraccion_mol_liquido[5],4)
        self.fraccion_mol_etilbenceno_Fase_vapor=round(self.fraccion_mol_vapor[5],4)
        
        self.fraccion_mol_dietilbenceno_Fase_liquida=round(self.fraccion_mol_liquido[6],4)
        self.fraccion_mol_dietilbenceno_Fase_vapor=round(self.fraccion_mol_vapor[6],4)
        

    def densidad(self):  ### SE DEFINEN LOS POLINOMIOS PARA EL CALCULO DE LA DENSIDAD FASE LIQUIDA (SACADOS DE ASPEN PLUS)
        
        if self.TC<=60:
            self.densidad_etileno_L=-2.37e-3*self.TC**3 + 0.289*self.TC**2-10.9*self.TC+ 339
        else:
            self.densidad_etileno_L=212.617
            
        if self.TC<=70:
            self.densidad_etano_L=1.51e-3*self.TC**3-0.116*self.TC**2 - 1.97*self.TC+ 402
        elif self.TC>70:
            self.densidad_etano_L=205.819
            
        if self.TC<=120:
            self.densidad_propileno_L=1e-5*self.TC**4 - 0.0023*self.TC**3 + 0.136*self.TC**2 - 4.13*self.TC + 552.37
        elif self.TC>120:
            self.densidad_propileno_L=225.72
            
        if self.TC<=280:
            self.densidad_benceno_L=-1.69e-5*self.TC**3 + 3.89e-3*self.TC**2-1.33*self.TC+ 903
        elif self.TC>280:
            self.densidad_benceno_L=9.27e-7*self.TC**4 - 1.52e-3*self.TC**3 + 0.922*self.TC**2 - 247*self.TC + 24878
            
        if self.TC<=300:
            self.densidad_tolueno_L=-2.22e-3*self.TC**2 - 0.594*self.TC + 878
        elif 300>=self.TC>400:
            self.densidad_tolueno_L=-7.97e-3*self.TC**3 + 0.882*self.TC**2 - 325*self.TC + 40163
        else:
            self.densidad_tolueno_L=298.53
            
        if self.TC<=350:
            self.densidad_etilbenceno_L=-2e-5*self.TC**3 + 0.0059*self.TC**2 - 1.4865*self.TC + 897.27
        elif self.TC>350:
            self.densidad_etilbenceno_L=4.18e-7*self.TC**4 - 7.42e-4*self.TC**3  + 0.493*self.TC**2 - 145*self.TC + 16277
            
        if self.TC<=390:
            self.densidad_dietilbenceno_L=-1.17e-5*self.TC**3 + 4.4e-3*self.TC**2 - 1.33*self.TC + 890
        elif self.TC>390:
            self.densidad_dietilbenceno_L=-6.5e-5*self.TC**3 + 0.0905*self.TC**2 - 41.9*self.TC + 6231
        
        #### DEFINE LAS DENSIDADES DE LOS COMPONENTES EN LA FASE VAPOR. SACADOS DE ASPEN PLUS
        
        self.densidad_etileno=3.69 -0.004747*self.P - 0.08858*self.T + 0.0001173*self.P**2 + 0.0008267*self.P*self.T + 0.0006316*self.T**2 -(2.993e-6)*self.P**2*self.T -(4.177e-6)*self.P*self.T**2 -(4.157e-6)*self.T**3 + 4.22e-4*self.P**2*self.T**2 -(2.056e-8)*self.P*self.T**3 + (3.854e-8)*self.T**4 -(6.723e-11)*self.P**2*self.T**3 + (1.232e-10)*self.P*self.T**4 -(1.27e-10)*self.T**5
        self.densidad_etano=0.002238 + 0.01233*self.P - (5.653e-5)*self.T + (1.125e-6)*self.P**2 - (4.515e-5)*self.P*self.T + (2.177e-6)*self.T**2 - (1.485e-8)*self.P**2*self.T + (1.696e-7)*self.P*self.T**2 - (4.226e-8)*self.T**3 + (8.265e-11)*self.P**2*self.T**2 - (5.57e-10)*self.P*self.T**3 + (3.07e-10)*self.T**4 - (1.176e-13)*self.P**2*self.T**3 + (9.471e-13)*self.P*self.T**4 - (7.299e-13)*self.T**5
        self.densidad_propileno=0.004475 + 0.01318*self.P - (7.861e-5)*self.T + (1.618e-06)*self.P**2 - (4.814e5)*self.P*self.T + (2.239e-6)*self.T**2 - (2.154e-8)*self.P**2*self.T + (1.858e-7)*self.P*self.T**2 - (4.7e-8)*self.T**3 + (1.214e-10)*self.P**2*self.T**2 - (6.518e-10)*self.P*self.T**3 + (3.643e-10)*self.T**4 - (2.542e-13)*self.P**2*self.T**3 + (1.177e-12)*self.P*self.T**4 - (8.986e-13)*self.T**5
        self.densidad_benceno=0.02613 + 0.0182*self.P - 0.0003562*self.T + (4.588e-6)*self.P**2 - (6.409e-5)*self.P*self.T + (3.292e-6)*self.T**2 - (6.466e-8)*self.P**2*self.T + (2.726e-7)*self.P*self.T**2 - (7.711e-8)*self.T**3 + (3.826e-10)*self.P**2*self.T**2 - (1.229e-9)*self.P*self.T**3 + (7.293e-10)*self.T**4 - (8.272e-3)*self.P**2*self.T**3 + (2.676e-12)*self.P*self.T**4 - (2.012e-12)*self.T**5
        self.densidad_tolueno=-2.102 + 0.04409*self.P + 0.1878*self.T + (9.719e-5)*self.P**2 -(0.0001563)*self.P*self.T - (0.002506)*self.T**2 + (1.608e-6)*self.P**2*self.T + (1.332e-5)*self.P*self.T**2 + (2.126e-5)*self.T**3 - (3.378e-8)*self.P**2*self.T**2 + (1.55e-8)*self.P*self.T**3 - (1.493e-7)*self.T**4 + (1.214e-10)*self.P**2*self.T**3 - (2.425e-10)*self.P*self.T**4 + (4.288e-10)*self.T**5
        self.densidad_etilbenceno=-18.93 + 0.2242*self.P + 0.5193*self.T - 0.0002407*self.P**2 - 0.006094*self.P*self.T - 0.002184*self.T**2 + (1.323e-5)*self.P**2*self.T + (3.735e-5)*self.P*self.T**2 - (6.819e-6)*self.T**3 - (1.364e-7)*self.P**2*self.T**2 + (6.459e-8)*self.P*self.T**3 - (2.512e-8)*self.T**4 + (3.891e-10)*self.P**2*self.T**3 - (6.403e-10)*self.P*self.T**4 + (3.186e-10)*self.T**5
        self.densidad_dietilbenceno=-9.936 + 0.2101*self.P - 0.2312*self.T - 0.0002825*self.P**2 + 0.001456*self.P*self.T + 0.001658*self.T**2 - (1.95e-6)*self.P**2*self.T - (3.087e-5)*self.P*self.T**2 + (1.777e-5)*self.T**3 + (9.056e-8)*self.P**2*self.T**2 - (4.353e-8)*self.P*self.T**3 - (4.894e-8)*self.T**4 - (3.793e-10)*self.P**2*self.T**3 + (7.399e-10)*self.P*self.T**4 -(2.603e-10)*self.T**5
            
        self.densidad_mezcla_vapor=self.densidad_etileno*self.fraccion_mol_etileno_Fase_vapor+self.densidad_etano*self.fraccion_mol_etano_Fase_vapor+self.densidad_benceno*self.fraccion_mol_benceno_Fase_vapor+self.densidad_tolueno*self.fraccion_mol_tolueno_Fase_vapor+self.densidad_etilbenceno*self.fraccion_mol_etilbenceno_Fase_vapor+self.densidad_dietilbenceno*self.fraccion_mol_dietilbenceno_Fase_vapor
        self.densidad_mezcla_liquida=self.densidad_etileno_L*self.fraccion_mol_etileno_Fase_liquida+self.densidad_etano_L*self.fraccion_mol_etano_Fase_liquida+self.densidad_benceno_L*self.fraccion_mol_benceno_Fase_liquida+self.densidad_tolueno_L*self.fraccion_mol_tolueno_Fase_liquida+self.densidad_etilbenceno_L*self.fraccion_mol_etilbenceno_Fase_liquida+self.densidad_dietilbenceno_L*self.fraccion_mol_dietilbenceno_Fase_liquida
        self.densidad_total_corriente=((self.Flujo_liquido*self.densidad_mezcla_liquida)/(self.corriente_alimento))+((self.Flujo_vapor*self.densidad_mezcla_vapor)/(self.corriente_alimento))
        
            
    def viscocidad(self): ### SE DEFINEN LOS POLINOMIO PARA EL CALCULO DE LA VISCOCIDAD (SACADOS DE ASPEN PLUS)
        
        
        ##### DEFINE LAS VISCOCIDADES PARA CADA COMPONENTE EN LA FASE VAPOR
        
        self.viscocidad_etileno_vapor=7e-15*self.T**3-2e-11*self.T**2+3e-8*self.T+9e-6       
        self.viscocidad_etano_vapor=6e-15*self.T**3-1e-11*self.T**2+3e-8*self.T+9e-6
        self.viscocidad_propileno_vapor=4e-15*self.T**3-1e-11*self.T**2+3e-8*self.T+8e-6
        self.viscocidad_benceno_vapor=5e-16*self.T**3-1e-12*self.T**2+3e-8*self.T+7e-6
        self.viscocidad_tolueno_vapor=3e-15*self.T**3+1e-11*self.T**2-2e-8*self.T+6e-6
        self.viscocidad_etilbenceno_vapor=3e-15*self.T**3-9e-12*self.T**2+2e-8*self.T+6e-6
        self.viscocidad_dietilbenceno_vapor=2e-15*self.T**3-8e-12*self.T**2+2e-8*self.T+5e-6        
        
        ##### DEFINE LAS VISCOCIDADES PARA CADA COMPONENTE EN LA FASE LIQUIDA
        
        self.viscocidad_etileno_liquido=-7e-13*self.T**3+7e-10*self.T**2-3e-7*self.T+5e-5
        self.viscocidad_etano_liquido=-2e-12*self.T**3+2e-9*self.T**2-5e-7*self.T+5e-5
        self.viscocidad_propileno_liquido=-3e-12*self.T**3+3e-9*self.T**2-9e-7*self.T+1e-4
        self.viscocidad_benceno_liquido=-2e-11*self.T**3+2e-8*self.T**2-6e-6*self.T+0.0008
        self.viscocidad_tolueno_liquido=-2e-11*self.T**3+2e-8*self.T**2-6e-6*self.T+0.0007
        self.viscocidad_etilbenceno_liquido=-2e-11*self.T**3+2e-8*self.T**2-6e-6*self.T+0.0008
        self.viscocidad_dietilbenceno_liquido=-4e-11*self.T**3+4e-8*self.T**2-1e-5*self.T+0.0013
        
        ### DEFINE LAS VISCOCIDADES DE LA MEZCLA TANTO VAPOR COMO LIQUIDA
        self.viscocidad_mezcla_vapor=self.viscocidad_etileno_vapor*self.fraccion_mol_etileno_Fase_vapor+self.viscocidad_etano_vapor*self.fraccion_mol_etano_Fase_vapor+self.viscocidad_propileno_vapor*self.fraccion_mol_propileno_Fase_vapor+self.viscocidad_benceno_vapor*self.fraccion_mol_benceno_Fase_vapor+self.viscocidad_etilbenceno_vapor*self.fraccion_mol_etilbenceno_Fase_vapor+self.viscocidad_dietilbenceno_vapor*self.fraccion_mol_dietilbenceno_Fase_vapor
        self.viscocidad_mezcla_liquido=self.viscocidad_etileno_liquido*self.fraccion_mol_etileno_Fase_liquida+self.viscocidad_etano_liquido*self.fraccion_mol_etano_Fase_liquida+self.viscocidad_propileno_liquido*self.fraccion_mol_propileno_Fase_liquida+self.viscocidad_benceno_liquido*self.fraccion_mol_benceno_Fase_liquida+self.viscocidad_tolueno_liquido*self.fraccion_mol_tolueno_Fase_liquida+self.viscocidad_etilbenceno_liquido*self.fraccion_mol_etilbenceno_Fase_liquida+self.viscocidad_dietilbenceno_liquido*self.fraccion_mol_dietilbenceno_Fase_liquida
        
        ##### DEFINE LA VISCOCIDAD PROMEDIO O TOTAL DE TODA LA CORRIENTE [VAPOR+LIQUIDA]
        self.viscocidad_total_corriente=((self.Flujo_liquido*self.viscocidad_mezcla_liquido)/(self.corriente_alimento))+((self.Flujo_vapor*self.viscocidad_mezcla_vapor)/(self.corriente_alimento))
        

    def entalpia(self): ### SE DEFINEN LOS POLINOMIO PARA EL CALCULO DE LA ENTALPIA (SACADOS DE ASPEN PLUS)
        
        ### DEFINE LA ENTALPIA PARA CADA COMPONENTE EN LA FASE LIQUIDA
        self.entalpia_etileno_liquido= 0.0393*(self.TC**2) + 50.828*self.TC + 46414
        self.entalpia_etano_liquido = 0.0004*(self.TC**3) - 0.2082*(self.TC**2) + 106.69*self.TC - 93068
        self.entalpia_propileno_liquido = 7e-11*(self.TC**6) - 1e-07*(self.TC**5) + 6e-05*(self.TC**4) - 0.0136*(self.TC**3) + 1.0856*(self.TC**2) + 151.83*self.TC + 1076.7
        self.entalpia_benceno_liquido = -5E-11*(self.TC**6) + 1E-07*(self.TC**5) - 8E-05*(self.TC**4) + 0.0249*(self.TC*3) - 3.129*(self.TC**2) + 263.78*self.TC + 45964
        self.entalpia_tolueno_liquido = 9E-11*(self.TC**6) - 1E-07*(self.TC**5) + 4E-05*(self.TC**4) - 0.0046*(self.TC**3) + 0.0528*(self.TC**2) + 167.21*self.TC + 9003.4
        self.entalpia_etilbenceno_liquido = 2E-10*(self.TC**6) - 3E-07*(self.TC**5) + 0.0001*(self.TC**4) - 0.0269*(self.TC*3) + 2.7035*(self.TC**2) + 79.495*self.TC - 14869
        self.entalpia_dietilbenceno_liquido = 2E-10*(self.TC**6) - 3E-07*(self.TC**5) + 0.0002*(self.TC**4) - 0.0416*(self.TC**3) + 4.8926*(self.TC**2) + 25.365*self.TC - 75753

        ### DEFINE LA ENTALPIA PARA CADA COMPONENTE EN LA FASE VAPOR

        self.entalpia_etileno_vapor= 0.0426*(self.TC**2) + 42.636*self.TC + 51328
        self.entalpia_etano_vapor= 0.0573*(self.TC**2) + 51.431*self.TC - 85254
        self.entalpia_propileno_vapor= 0.0652*(self.TC**2) + 64.204*self.TC + 18434
        self.entalpia_benceno_vapor= 0.111*(self.TC**2) + 86.483*self.TC + 80149
        self.entalpia_tolueno_vapor= -7E-05*(self.TC**3) + 0.1858*(self.TC**2) + 98.424*self.TC + 47161
        self.entalpia_etilbenceno_vapor= -1E-04*(self.TC**3) + 0.2269*(self.TC**2) + 119.39*self.TC + 26225
        self.entalpia_dietilbenceno_vapor= -1E-04*(self.TC**3) + 0.2735*(self.TC**2) + 167.15*self.TC - 27343
        
        ### DEFINE LA ENTALPIA DE LA MEZCLA TANTO DE VAPOR COMO LIQUIDA
        self.entalpia_mezcla_vapor=self.entalpia_etileno_vapor*self.fraccion_mol_etileno_Fase_vapor+self.entalpia_etano_vapor*self.fraccion_mol_etano_Fase_vapor+self.entalpia_propileno_vapor*self.fraccion_mol_propileno_Fase_vapor+self.entalpia_benceno_vapor*self.fraccion_mol_benceno_Fase_vapor+self.entalpia_tolueno_vapor*self.fraccion_mol_tolueno_Fase_vapor+self.entalpia_etilbenceno_vapor*self.fraccion_mol_etilbenceno_Fase_vapor+self.entalpia_dietilbenceno_vapor*self.fraccion_mol_dietilbenceno_Fase_vapor
        self.entalpia_mezcla_liquido=self.entalpia_etileno_liquido*self.fraccion_mol_etileno_Fase_liquida+self.entalpia_etano_liquido*self.fraccion_mol_etano_Fase_liquida+self.entalpia_propileno_liquido*self.fraccion_mol_propileno_Fase_liquida+self.entalpia_benceno_liquido*self.fraccion_mol_benceno_Fase_liquida+self.entalpia_tolueno_liquido*self.fraccion_mol_tolueno_Fase_liquida+self.entalpia_etilbenceno_liquido*self.fraccion_mol_etilbenceno_Fase_liquida+self.entalpia_dietilbenceno_liquido*self.fraccion_mol_dietilbenceno_Fase_liquida
        
        ### DEFINE LA ENTALPIA PROMEDIO O TOTAL DE TODA LA CORRIENTE [VAPOR+LIQUIDA]
        self.entalpia_total_corriente=((self.Flujo_liquido*self.entalpia_mezcla_liquido)/(self.corriente_alimento))+((self.Flujo_vapor*self.entalpia_mezcla_vapor)/(self.corriente_alimento))
        
        
    def imprimir_en_pantalla(self): #### MODULO PARA VER EN PANTALLA LAS PROPIEDADES TERMODINAMICAS DE LA MEZCLA

        print("-----------------------------------------------------------------------------")
        print("Las raices Z son : ", self.raices)
        print("-----------------------------------------------------------------------------")        
        print("El factor de compresibilidad de la mezcla [vapor] es : ",self.Z_vapor)
        print("-----------------------------------------------------------------------------")
        print("El factor de compresibilidad de la mezcla [liquida] es : ",self.Z_liquido)
        print("-----------------------------------------------------------------------------")
        print("Los coeficientes K para cada componente son : ",self.coeficientes_K)
        
        
        print("______________________________________________________________________________")
        print("PARAMETROS CALCULADOS [RACHFORD RICE]")
        print("______________________________________________________________________________")
        print("El flujo molar de la corriente de alimento [F] es : ",self.corriente_alimento,"[Kmol/h]")
        print("-----------------------------------------------------------------------------")
        print("El flujo molar de la fase liquida [L] es : ", self.Flujo_liquido,"[Kmol/h]")
        print("-----------------------------------------------------------------------------")
        print("El flujo molar de la fase vapor [V] es : ", self.Flujo_vapor,"[Kmol/h]")
        print("-----------------------------------------------------------------------------")
        print("Las fracciones molares para el etileno en la fase vapor-liquida  son : " ,self.fraccion_mol_etileno_Fase_vapor,"-",self.fraccion_mol_etileno_Fase_liquida)
        print("-----------------------------------------------------------------------------")
        print("Las fracciones molares para el etano en la fase vapor-liquida son : " ,self.fraccion_mol_etano_Fase_vapor,"-",self.fraccion_mol_etano_Fase_liquida)
        print("-----------------------------------------------------------------------------")
        print("Las fracciones molares para el propileno en la fase vapor-liquida son : " ,self.fraccion_mol_propileno_Fase_vapor,"-",self.fraccion_mol_propileno_Fase_liquida)
        print("-----------------------------------------------------------------------------")
        print("Las fracciones molares para el benceno en la fase vapor-liquida son : " ,self.fraccion_mol_benceno_Fase_vapor,"-",self.fraccion_mol_benceno_Fase_liquida)
        print("-----------------------------------------------------------------------------")
        print("Las fracciones molares para el tolueno en la fase vapor-liquida son : " ,self.fraccion_mol_tolueno_Fase_vapor,"-",self.fraccion_mol_tolueno_Fase_liquida)
        print("-----------------------------------------------------------------------------")
        print("Las fracciones molares para el etilbenceno en la fase vapor-liquida son : " ,self.fraccion_mol_etilbenceno_Fase_vapor,"-",self.fraccion_mol_etilbenceno_Fase_liquida)
        print("-----------------------------------------------------------------------------")
        print("Las fracciones molares para el dietilbenceno en la fase vapor-liquida son : " ,self.fraccion_mol_dietilbenceno_Fase_vapor,"-",self.fraccion_mol_dietilbenceno_Fase_liquida)
        print("______________________________________________________________________________")
        
        
        print("*****************DENSIDAD****************")
        print("______________________________________________________________________________")
        print("-----------------------------------------------------------------------------")
        print("La densidad de la mezcla de gases en la fase vapor es : ", self.densidad_mezcla_vapor,"[Kg/m3]")
        print("-----------------------------------------------------------------------------")
        print("La densidad de la mezcla de gases en la fase liquida es : ", self.densidad_mezcla_liquida,"[Kg/m3]") 
        print("-----------------------------------------------------------------------------")
        print("La densidad total de la corriente es : ", self.densidad_total_corriente,"[Kg/m3]")
        print("______________________________________________________________________________")
        

        print("************VISCOCIDAD********************")
        print("-----------------------------------------------------------------------------")
        print("La viscocidad de la mezcla  en la fase vapor es : ", self.viscocidad_mezcla_vapor,"[Kg/m s]")
        print("-----------------------------------------------------------------------------")
        print("La viscocidad de la mezcla  en la fase liquida es : ", self.viscocidad_mezcla_liquido,"[Kg/m s]") 
        print("-----------------------------------------------------------------------------")
        print("La viscocidad total de la corriente es : ", self.viscocidad_total_corriente,"[Kg/m s]")
        
        
        print("______________________________________________________________________________")
        print("****************ENTALPIA*********************")
        print("______________________________________________________________________________")
        print("La entalpia de la corriente fase vapor es : ",self.entalpia_mezcla_vapor,"[KJ/Kmol]")
        print("-----------------------------------------------------------------------------")
        print("La entalpia de la corriente fase liquida es : ",self.entalpia_mezcla_liquido,"[KJ/Kmol]") 
        print("-----------------------------------------------------------------------------")
        print("La entalpia total de la corriente es : ",self.entalpia_total_corriente,"[KJ/Kmol]")
        
                
#### MODULOS REALIZADOS [NO BORRAR NINGUNO]
        
        
flash=Simulador_etilbenceno()
flash.parametros()
flash.parametros_Ecuacion_Estado()
flash.factor_compresibilidad()
flash.raices_reales_imaginarias()
flash.coeficiente_fugacidad()
flash.Rachford_Rice()
flash.fracciones_salida_Flash()
flash.densidad()
flash.viscocidad() 
flash.entalpia() 
flash.imprimir_en_pantalla() 

 ### LINEA DE CODIGO PARA CORRER MODO ONLINE (ideone).online 
sys.stdin.readline()