%Programa demostrativo de c�lculo de altura de un equipo con transferencia de masa
%que est� operando en:

%                           CO-CORRIENTE

% NOTA: si la altura de la torre es insuficiente, el programa corre pero al
%final muestra en pantalla este error:

%  Undefined function or variable 'XA_OK'.

%  Error in SimulacionparaAlturaEquipoTdeM_Cocorriente (line 207)
%        plot(XAentrada,YAentrada,'dk','LineWidth',2),plot(XA_OK,YA_OK,'dk','LineWidth',2) 

% No se debe alarmar por eso. Significa que hace falta altura y como no se
% encontr� la limpieza deseada del Refinado, el programa NO puede marcar el
% punto final de la LdeO. Aumente longitud de torre para ver el resultado OK.

%Se recalca que el modelo es en Estado Estacionario para todas la torre y para
%sus particiones. Es decir, se asume que las concentraciones y condiciones
%de cada rodaja o partici�n son las que alcanza la torre cuando lleva
%suficiente tiempo de estar operando.

%Se asume que el Refinado R es la fase de inter�s, en este caso un Gas, que esta en eje y.
%La fase de servicio es el Extracto E es un l�quido, que est� en el eje x.
%   R=GAS     E=L�quido

%Se tomar� como si el equipo fuera una Torre Empacada con anillos Raschig (vertical) de altura
%final z a calcular para lograr el nivel de separaci�n dado para fase R. A esta torre las dos
%corrientes entran por la parte superior (gas y l�quido) y salen por el fondo.

%Todo se asume en unidades MOLARES

clear all
close all
clc

%Longitud total a simular de la altura de la torre. 
%Dada en metros. Se sugiere entrar al principio valores con mucha holgura 
%para ver estado estacionario del equipo (por ejemplo 5 metros).
LongitudMaxima=input('Entre la altura de la torre en metros  = ');

%Datos del equipo
%Concentraciones en el Refinado: entrada=Dato  salida=Esperada
%Primero medidas en laboratorio como fracci�n respecto de la corriente como
%un todo. Recuerden Formulaciones de conversi�n:   XA=xA/(1-xA) o  xA=XA/(1+XA)

yAentrada=0.4;      %Fracci�n molar respecto de R como un todo Gas Eje y
yAsalida=0.22;      %Fracci�n molar respecto de R como un todo Gas Eje y
%Conversi�n de concentraciones a fracci�n molar soluto respecto ste puro
YAentrada=yAentrada/(1-yAentrada);      %Fracci�n molar del soluto respecto del solvente libre de soluto
YAsalida=yAsalida/(1-yAsalida);        %Fracci�n molar del soluto respecto del solvente libre de soluto

%Concentraciones en el Extracto: entrada=Dato  salida=Calculada, por eso NO
%se conoce para declararla como un dato.
%Primero medidas en laboratorio como fracci�n respecto de la corriente como
%un todo
xAentrada=0.05;        %Fracci�n molar respecto de E como un todo
%Conversi�n de concentraciones a fracci�n molar soluto respecto ste puro
XAentrada=xAentrada/(1-xAentrada);      %Fracci�n molar del soluto respecto del solvente libre de soluto

%Muestra al usuario valores en unidades de concentraci�n respecto del
%solvente libre de soluto
fprintf('\n\nCONCENTRACIONES EN FRACCI�N RESPECTO DEL SOLVENTE PURO.\n\n');
Fase = 'Refinado';
Punto = ' Entrada';
Valor=num2str(YAentrada);
fprintf('A la %s el %s contiene %s kgA/kgSteR_Puro.\n', Punto, Fase, Valor);

Punto = ' Salida';
Valor=num2str(YAsalida);
fprintf('A la %s el %s contiene %s kgA/kgSteR_Puro.\n', Punto, Fase, Valor);

Fase = 'Extracto';
Punto = ' Entrada';
Valor=num2str(XAentrada);
fprintf('A la %s el %s contiene %s kgA/kgSteR_Puro.\n\n', Punto, Fase, Valor);
fprintf('PULSE ENTER PARA CONTINUAR.\n\n');
pause

%Datos de las corrientes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Refinado=GAS
%Flux molar de entrada de la fase Refinado
fluxMol_REntrada=0.05;      %kmolMixR/m2/s
%Flux del Solvente puro en la fase Refinado
fluxRsolvente=fluxMol_REntrada*(1-YAentrada); %Es constante
%Masa molecular de la fase refinado
MMR=11;         %kg/kmol
%Viscosidad de la fase Refinado a las condiciones de la entrada
muREntrada=1e-5;    %kg/m-s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Extracto=L�quido
%Flux molar de entrada de la fase Extracto
fluxMol_EEntrada=0.065;      %kmol/m2-s
%Flux del Solvente puro en la fase Extracto
fluxEsolvente=fluxMol_EEntrada*(1-XAentrada); %Es constante
%Masa molecular de la fase Extracto
MME=260;         %kg/kmol
%Viscosidad de la fase Extracto
muEEntrada=2e-3;    %kg/m-s

%Pendiente de la L�nea de Operaci�n LdeO calculada con los fluxes
mLdeO=-fluxEsolvente/fluxRsolvente;


%Datos del Equipo
Poros=0.75;     %Porosidad para este tipo de empaque
CorrecHoldUpLiq=0.0333; %Correcci�n de porosidad por retenido de l�quido
%Porosidad Operativa de la torre (descontando ocupaci�n del agua)
PorosOper=Poros-CorrecHoldUpLiq;
%D�ametro del empaque usado en la torre
ds=0.0472;      %m
%�rea para la transferencia de masa dada en m2/m3 de empaque en la torre
aM=37.4;        %m2/m3
%Recu�rdese que NO se requiere el di�metro de la torre puesto que todo se
%calcula con los fluxes, por lo que el di�metro se halla con el flujo molar
%de una de las fases como dato.

   
%Altura de torre que se tomar� en cada partici�n
PasoEspacial=0.05;       %metros. Verificar que sea mayor que ds di�metro del empaque, en este caso >0.0472m

%Numero de pasos del ciclo for para reconstruir el tiempo total. Notese que
%se redondea hacia arriba (operador ceil=cielo) para no tener numero de pasos fraccionario.
NumPasos=ceil(LongitudMaxima/PasoEspacial);


%Inicializa variables de entrada de las dos fases en el vector que les corresponde
%estar lista en la primera iteraci�n
ConcR(1)=YAentrada;
ConcE(1)=XAentrada;
fluxMol_R(1)=fluxMol_REntrada;
fluxMol_E(1)=fluxMol_EEntrada;

%Ciclo de recorrido de la torre desde la entrada hasta una altura suficientemente grande para poder ver
%la evoluci�n de las variables de concentraci�n hasta su estado estacionario.
for i=1:NumPasos    
   %Calcular coeficiente local englobante de TdeM, que se toma como el del gas asumiendo que el coeficiente 
   %local para el l�quido es muy grande por lo que el limitante es el de la fase R=gas, pero podr�a ser 
   %al rev�s.
   %Prec�lculo de variables actualizadas para el coeficiente. Actualiza la Viscosidad de
   %la fase R de acuerdo con su concentraci�n de soluto, pues afecta el coeficiente de TdeM
   muR=muREntrada+(ConcR(i)/100)^2;
   
   %Actualiza coeficiente para la secci�n tratada en unidades de kmol/m2-s
   %pues la concentraciones son fracciones adimensionales
   KY(i)=fluxMol_R(i)*0.15*((ds*(fluxMol_R(i)*MMR))/(muR*(1-PorosOper)))^-0.36;
   
   %Calcula la concentraci�n ficticia correspondiente al punto actual la cual se despeja de
   %la CDE: Y*=f(X). N�tese que se debe conocer la concentraci�n en la fase Extracto X=ConcE(i) 
   %como dato para entrar a la CDE y hallar la concentraci�n ficticia
   ConcRFicticia(i)= 9.524*(ConcE(i))^3-7.302*(ConcE(i))^2 + 2.064*(ConcE(i));  
   %Calcula el Flux de transferencia de masa NA usando concentraciones
   %ficticias (*)
   NA(i)=KY(i)*(ConcR(i)-ConcRFicticia(i));   %kmol/m2-s por m2 de �rea de Transferencia de Masa
   %Calcula la cantidad transferida en ese tramo de torre
   NAtramo(i)=NA(i)*aM*PasoEspacial;        %kmol/m2-s * m2/m3  [=] kmol/m2-s por m2 de �rea de flujo
   
      %Actualiza valores de propiedades de las fases para el inicio del tramo siguiente
   %Para la fase de inter�s: Refinado
   %Aplica transferencia de masa al flux, pues NA del tramo tambi�n es flux
   fluxMol_R(i+1)=fluxMol_R(i)-NAtramo(i);  %flux en kmol/m2-s
  
   %Actualiza concentraci�n, usando la expresi�n de balance de masa por
   %componente del soluto en la secci�n tratada:
   %                                            Y(i)*Rsolvente=Y(i+1)*Rsolvente+TdeMasa
   ConcR(i+1)=ConcR(i)-(NAtramo(i)/(fluxRsolvente));
   
   %Actualiza el Flux Molar en este caso para la fase Extracto, que recibe soluto
   fluxMol_E(i+1)=fluxMol_E(i)+NAtramo(i);
   
   %Actualiza concentraci�n, usando la expresi�n de balance de masa por
   %componente del soluto en la secci�n tratada:
   %X(i)*Esolvente+TdeMasa=X(i+1)*Esolvente
   ConcE(i+1)=ConcE(i)+(NAtramo(i)/(fluxEsolvente));
   
   %Ubica la concentraci�n del Extracto que corresponde a una concentraci�n
   %del Refinado aproximadamente igual a la deseada de salida
   if ConcR(i+1)>0.995*YAsalida & ConcR(i+1)<1.005*YAsalida
       %Valor de concentraci�n del Extracto cuando el Refinado esta listo
       XA_OK=ConcE(i+1);
       YA_OK=ConcR(1)+mLdeO*(XA_OK-ConcE(1));
   end
   
   %Crea vector de concentraci�n deseada de salida del Refinado para graficarla
   ConcRDeseada(i)=YAsalida;
   
   %Crea vector de eje altura de la torre para ponerlo en el eje x del gr�fico
   AlturaTorre(i)=i*PasoEspacial;
      
end

%Ajusta longitud de todos los vectores a la longitud total de torre probada. Solo lo hace para efecto gr�ficos
%Alarga vectores que se quedaron en indice i adicion�ndoles un valor
ConcRFicticia(i+1)= 9.524*(ConcE(i+1))^3-7.302*(ConcE(i+1))^2 + 2.064*(ConcE(i+1));
ConcRDeseada(i+1)=YAsalida;
AlturaTorre(i+1)=(i+1)*PasoEspacial;


figure(1), plot(AlturaTorre,ConcR,'*r'), grid on, hold on, plot(AlturaTorre,ConcR,'r'), plot(AlturaTorre,ConcE,'*b'),...
                plot(AlturaTorre,ConcRFicticia,'*g'), plot(AlturaTorre,ConcRDeseada,'.k'),...
                title('Concs. Refinado: Rojo, Extracto: Azul, Fictica: Verde y Deseada R: Negro')
pause

%Grafica la CDE y la LdeO
for j=1:300
    %CDE
    ConcEFig(j)=j*0.001;
    CDE(j)=9.524*(ConcEFig(j))^3-7.302*(ConcEFig(j))^2 + 2.064*(ConcEFig(j));
end


%Grafica la CDE y la LdeO

for j=1:120
    ConcELdeO(j)=ConcE(1)+(j-1)*0.001;
    Y_LdeO(j)=ConcR(1)+mLdeO*(ConcELdeO(j)-ConcE(1));
end

figure(2), plot(ConcEFig,CDE,'.k'), grid on, hold on, plot(ConcEFig,CDE,'k'), plot(ConcELdeO,Y_LdeO,'.g'), plot(ConcELdeO,Y_LdeO,'g'),...
        plot(ConcE,ConcR,'.r'), title('CDE: Negro, LdeO: Verde y Evoluci�n del proceso: Rojo'),...
        plot(XAentrada,YAentrada,'dk','LineWidth',2),plot(XA_OK,YA_OK,'dk','LineWidth',2)