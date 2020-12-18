%%%%%%%% IDENTIFICACACION DEL PROCESO
% Este programa permite hallar el estado estacionario de un equipo de
% deshumidificacion
% operando en Contra-Corriente. Dicho estado estacionario se puede usar solo en este proceso.

% Autores: JULIANA- DAVID-JOSE
% �ltima actualizaci�n: Abri, 2018

% Nomenclatura usada y esquem�tico de corrientes:
% R: Refinado. SdeP I. Solucion de aire con benceno. R soluci�n
% completa, RS solvente en la soluci�n de la fase R
% E: Extracto. SdeP II. Benceno liquido(GOTAS). E Benceno liquido
% X_A: Concentracion de benceno en el Extracto [kmol benceno/kmol
% benceno]que va hacer igual a uno  por tener benceno como soluto y
% reefinado
% Y_A: Concentraci�n de A la fase Refinado [kmol benceno /kmol aire]

%    X_A(0)      Y_A(1)    
%       -----------    ---- Plano 1, partici�n 1:  Entrada del Extracto (Salida de R)
%       |         |
%       |         |
%       |  Torre  |
%       |         |
%       |         |
%       -----------    ---- Plano 2, partici�n NP:  Entrada del Refinado (Salida de E)
%    X_A(NP)   Y_A(NP+1) 

% En t�rminos de vectores, la posici�n j siempre indica flujos y concentraciones correspondientes
% a la partici�n j, pero en sus salidas, NO en las entradas. Por eso, al graficar X(1),Y(1), 
% no se ver� la entrada del Refinado a la torre sino la salida de R de la primera partici�n
% con la salida de la primera partici�n o salida de E de la torre. Debe prestarse atenci�n
% a esto porque estas parejas forman la L�nea de Operaci�n, pero no son las que se usan para el
% c�lculo de las propiedades dentro de la partici�n, pues una es salida y la otra es entrada.

%%%%%%%%% AJUSTE DE LA PLATAFORMA DE PROGRAMACI�N
clear all 
close all
clc
format long

%%%%%%%% INICIALIZACI�N DE VARIABLES GENERALES
% Altura de la partici�n a tomar en la torre
dL = 0.05;                          % [m] 
% Altura �til de la torre
L_T = 3.1;                   %[m]
% N�mero de particiones hechas en la torre para resolver el modelo
NP = ceil(L_T/dL);       % [Adim] 
% Di�metro de la torre
D_T = 1;                        %[m]
% Tolerancia del error en la iteraci�n para convergencia
Tol_error = 1e-8;%
% Valor m�ximo de ciclos de iteraci�n de convergencia
MaxContador = 100;

%%%%%%%% DECLARACI�N DE CONSTANTES DEL MODELO
% Aceleraci�n de la gravedad
g = 9.8;                        %[m/s2]
% Tensi�n superficial de la fase dispersa: soluci�n de Tolueno
sigma_T = 0.0291;                %[N/m]
% Viscosidad del Tolueno puro
mu_T = 5.24825e-4;               %[Pa-s]
% Viscosidad del Agua pura
mu_W = 8e-4;                     %[Pa-s]
%Viscosidad de la Nicotina pura
mu_N = 3.442e-3;                 %[Pa-s]
% Masa molecular del Tolueno 
MM_T = 92.138;                   %[kg/kmol]
% Masa molecular del Agua
MM_W = 18;                       %[kg/kmol]
% Masa molecular de la Nicotina
MM_N = 162.232;                  %[kg/kmol]
% Difusividad de la Nicotina en Agua
Dif_N_W = 6.8753e-10;                  
% Difusividad de la Nicotina en Tolueno
Dif_N_T = 1.62033e-9; 
% Las siguientes propiedades dependen de la temperatura, por lo que en caso
% de operar a otra T se deber�n ajustar a los valores correctos.
% Densidad del Tolueno a 303K
rho_T = 857.7;                   %[kg/m3]
% Densidad del Agua a 303K
rho_W = 995.7;                   %[kg/m3]
% Densidad de la Nicotina a 303K
rho_N = 1010;                    %[kg/m3]

%%%%%%%%% INICIALIZACI�N DE VARIABLES ESCALARES DEL MODELO
% Di�metro de la gota de la fase dispersa: soluci�n de Nicotina en Tolueno,
% dada en mil�metros pero convertida a metros.
D_g = 1/1000;                      % [m]
% Flujo molar de solvente en el Refinado: Agua, que es cte en toda la torre
% dado en moles por segundo pero convertido a kmoles
dotn_RS = 0.2;                % [kmol/s] de Agua
% Flujo molar de se solvente en el Extracto: Tolueno, cte en toda la torre
% dado en moles por segundo pero convertido a kmoles
dotn_ES = 0.0158;            % [kmol/s] de Tolueno
% Concentraciones en el Plano 1, Superior:
% Concentraci�n de Entrada (dato) del Refinado. Corresponde a la posici�n Cero 0
X_A_Ent = 0.0087;                %[kmol Nicotina/kmol RS]
% Concentraci�n de Salida Esperada del Extracto
Y_A(1) = 0.1;                 %[kmol Nicotina/kmol ES]
% Concentraciones en el Plano 2, Inferior:
% Concentraci�n de Entrada (dato) del Extracto. Corresponde a la posici�n NP+1
Y_A_Ent = 0.005;              %[kmol Nicotina/kmol ES]
% Concentraci�n de Salida Esperada del Refinado
X_A(NP) = 0.0012;               %[kmol Nicotina/kmol RS]

% C�lculos usando las inicializaciones previas.
% �rea de flujo total a torre vac�a
A_FT = pi*(D_T^2)/4;            %[m2]
% Volumen de una partici�n
VP = A_FT*dL;                   %[m3]
% Intervalo total de concentraciones del Refinado (con la Conc. de salida Esperada)
Int_X_A = X_A_Ent - X_A(NP);    %[kmol Nicotina/kmol RS]
% Calcula la pendiente de la l�nea de operaci�n que se espera de la torre
m_LdeO = (Y_A_Ent - Y_A(1))/(X_A(NP) - X_A_Ent);    %[kmol RS/kmol ES]
% Volumen de una gota
V_g = (1/6)*pi*D_g^3;           %[m3]
% �rea superficial de una gota
As_g = pi*D_g^2;                %[m2]

%%%%%%%%% INICIALIZACI�N DE VARIABLES VECTORIALES DEL MODELO
% Genera puntos para graficar la l�nea de operaci�n LdeO. Asigna la primera
% pareja, correspondiente a las concentraciones en el Plano 1
X_LdeO(1)= X_A_Ent;
Y_LdeO(1) = Y_A(1);
% Ciclo for para generar el resto de puntos de la l�nea de operaci�n LdeO
for l = 1:NP
    % Valor de concentraci�n supuesta para el Refinado: con espaciado uniforme, a la entrada
    % de la partici�n l+1 que es lo mismo que la salida de la partici�n l.
    X_LdeO(l+1) = X_LdeO(l)-(Int_X_A/NP);
    % Calcula el valor correspondiente de la concentraci�n de Extracto, que es la de salida de
    % la partici�n l+1, usando la ecuaci�n de la l�nea de operaci�n
    Y_LdeO(l+1) = m_LdeO*(X_LdeO(l+1) - X_A(NP)) + Y_A_Ent;
    % N�tese que la LdeO contiene siempre la coordenada X de la entrada a la partici�n, 
    % mientras que la coordenada Y corresponde a la salida de esa partici�n. Es decir,
    % no son puntos de salidas de la partici�n tratada.
end

% Inicializa el vector de concentraciones finales a obtener del Refinado, en el cual
% la primera posici�n es conocida e igual a la concentraci�n de entrada de R dada
X_A(1) = X_A_Ent;

%%%%%%%%%% RESUELVE TODA LA TORRE POR PARTICIONES HASTA LOGRAR LA CONVERGENCIA
for j = 1:NP  % Ciclo for para la variaci�n de la Posici�n de la Partici�n en la Torre
    % Asigna una concentraci�n semilla para la iteraci�n en la salida del Refinado
    % de la partici�n actual, usando una repartici�n uniforme de distribuci�n de 
    % concentraci�n. Parte de la concentraci�n conocida de la entrada de R.
    X_Asup(j) = X_A_Ent - j*(Int_X_A/NP);
    
    if j==1
        % Asigna la concentraci�n de Refinado en la corriente de entrada a la partici�n
        % para evitar el problema del �ndice cero en Matlab para los vectores
        X_A_EntPart(j) = X_A_Ent; % La primera vez es la de entrada
    else
        % Asigna la concentraci�n de Refinado en la corriente de entrada a la partici�n
        % para el resto de particiones, que no presentan el problema del �ndice cero en Matlab
        X_A_EntPart(j) = X_A(j-1); % Las otras veces es la anterior
    end    
    
    % La concentraci�n de la fase Extracto a la salida de esta partici�n se conoce para
    % la primera partici�n, pues es la concentraci�n de salida esperada. Para las dem�s
    % particiones ser� un valor que la iteraci�n previa ya calcul� como el valor 
    % definitivo del estado estacionario.
    if j==1
        Y_A(j) = Y_A(1);
    end
          
    % Ciclo while de iteraci�n para buscar convergencia en la concentraci�n del 
    % Refinado a la salida de cada partici�n del equipo.
    % Inicializa el error en un valor dos veces la torelacia para entrar a la iteraci�n
    errorEst = 2*Tol_error;
    % Inicializa contador para el n�mero de iteraciones m�ximas
    contador = 0; 
    while errorEst > Tol_error && contador < MaxContador
        % Bloque de c�lculo del Flux de trasferencia de masa.
        % Como se requiere tener las concentraciones en las unidades adecuadas para
        % actualizar las propiedadess de las fases, hace las conversiones necesarias.
        % Conversi�n a fracci�n molar respecto de la fase como un todo
        x_Asup(j) = X_Asup(j)/(1+X_Asup(j));
        y_A(j) = Y_A(j)/(1+Y_A(j)); 
        
        % Conversi�n a fracci�n m�sica respecto del solvente puro
        Xw_Asup(j) = X_Asup(j)*MM_N/MM_W;
        Yw_A(j) = Y_A(j)*MM_N/MM_T;
    
        % Conversi�n a fracci�n m�sica respecto de la fase como un todo
        xw_Asup(j) = Xw_Asup(j)/(1+Xw_Asup(j));
        yw_A(j) = Yw_A(j)/(1+Yw_A(j));     
    
        % Primero se calculan todas las propiedades de las fases requeridas para
        % hallar los coeficientes. Todo se trata como mezcla ideal.
        % Refinado R, SdeP I
        % Densidad de la mezcla R
        rho_R(j) = 1/((xw_Asup(j)/rho_N)+((1-xw_Asup(j))/rho_W));   %[kg/m3]
        % Masa molecular de la mezcla R
        MM_R(j) = MM_W*(1-x_Asup(j)) + MM_N*x_Asup(j);              %[kg/kmol]
        % Viscosidad de la mezcla R
        mu_R(j) = exp(x_Asup(j)*log(mu_N)+(1-x_Asup(j))*log(mu_W)); %[Pa-s]
    
        % Extracto E, SdeP II
        % Densidad de la mezcla E
        rho_E(j) = 1/((yw_A(j)/rho_N)+((1-yw_A(j))/rho_T));         %[kg/m3]
        % Masa molecular de la mezcla E
        MM_E(j) = MM_T*(1-y_A(j)) + MM_N*y_A(j);                    %[kg/kmol]
        % Viscosidad de la mezcla E
        mu_E(j) = exp(y_A(j)*log(mu_N)+(1-y_A(j))*log(mu_T));       %[Pa-s]
        
        % Como va a requerir flujos volum�tricos de las fases para los
        % c�lculos, los halla con los datos de flujos molares existentes.
        % Calcula el flujo total de Refinado como la suma del solvente y el soluto 
        dotn_R(j) = dotn_RS + dotn_RS*X_Asup(j);        %[kmol/s]
        % Lo convierte a flujo m�sico con la masa molecular del Refinado en ese punto
        dotm_R(j) = dotn_R(j)*MM_R(j);                  %[kg/s]
        % Ahora lo conviere a flujo volum�trico usando la densidad actual del Refinado
        dotV_R(j) = dotm_R(j)/rho_R(j);                 %[m3/s]
        
        % Calcula el flujo total de Extracto como la suma del solvente y el soluto 
        dotn_E(j) = dotn_ES + dotn_ES*Y_A(j);           %[kmol/s]
        % Lo convierte a flujo m�sico con la masa molecular del Extracto en ese punto
        dotm_E(j) = dotn_E(j)*MM_E(j);                  %[kg/s]
        % Ahora lo conviere a flujo volum�trico usando la densidad actual del Refinado
        dotV_E(j) = dotm_E(j)/rho_E(j);                 %[m3/s]
        
        % Calcula el coeficiente local de TdeM para la Fase Refinado, con concentraciones X
        % C�lculo del n�mero de Reynolds de la gota, pero primero halla todo lo previo
        % Fracci�n volum�trica del Extracto o hold-up, fase dispersa, en la torre
        FV_E(j) = dotV_E(j)/(dotV_R(j) + dotV_E(j));
        % Velocidad de la fase Refinado en la torre
        v_R(j) = dotV_R(j)/(A_FT*(1-FV_E(j)));          %[m/s]
        % Velocidad de la fase Extracto en la torre
        v_E(j) = dotV_E(j)/(A_FT*FV_E(j));              %[m/s]
        % Velocidad de deslizamiento entre las fases.
        v_s(j) = v_R(j) + v_E(j);                       %[m/s]
        % N�mero de Reynolds de la gota, porque as� lo pide la correlaci�n
        Re_g(j) = ((D_g*v_s(j)*rho_R(j))/mu_R(j));
        
        % C�lculo del n�mero de Schmidt
        Sc_R(j) = mu_R(j)/(rho_R(j)*Dif_N_W);  
        
        % Coeficiente propio de la correlaci�n de Skellan
        K_Skll = (Re_g(j)^(1/8))*((mu_R(j)/mu_E(j))^(1/4))*((mu_R(j)*v_s(j)/sigma_T)^(1/6));
        % Factor de correcci�n por forma
        F_Skll = 0.281+ 1.615*K_Skll+3.73*K_Skll^2-1.874*K_Skll^3;
        % Coeficiente local de TdeM para la fase continua, el Refinado
        k_X(j) = (Dif_N_T/D_g)*(2+0.463*Re_g(j)^(0.484)*Sc_R(j)^(0.339)*...
                    ((D_g*g^(1/3))/(Dif_N_W^(2/3)))^(0.072))*F_Skll;
        
        % Calcula coeficiente local de TdeM  para Fase Extracto, con concentraciones Y
        k_Y(j) = 0.00375*v_s(j)/(1+(mu_E(j)/mu_R(j)));
        
        % Calcula la concentraci�n en la interfase, para lo cual primero halla la
        % pendiente de la l�nea de Fuerza Impulsora (FI)
        m_FI = -k_Y(j)/k_X(j);
        % Resuelve la ecuaci�n resultante del cruce de la FI y la CDE usando el punto actual
        % Usa la funci�n que encuentra el cero o la ra�z de una funci�n no lineal
        % Declara la funci�n Cruce de FI y CDE en un .m externo y declara una funci�n
        % muda FunConArg que llama a la funci�n original para pasarle los par�metros
        FunConArg=@(x)CruceFI_CDE(x,X_Asup(j),Y_A(j),m_FI);
        % Da un punto inicial para buscar la soluci�n
        if j==1
            % Para la primera partici�n da un valor arbitrario que est� cerca del
            % valor esperado. Este valor se determinan mirando la LdeO y la CDE.
            x_0 = 0.001;
        else
            % Para las dem�s particiones toma como semilla el valor anterior
            % buscando mantener la continuidad de la salida 
            x_0 = X_Aint(j-1);
        end
        % Llama la funci�n que encuntra la ra�z o el cero de la funci�n no lineal
        X_Aint(j) = fzero(FunConArg,x_0);
        
        % Calcula el flux de TdeM usando la ecuaci�n con coeficiente local y las
        % concentraciones asociadas con la fase R
        N_A(j) = k_X(j)*(X_Asup(j) - X_Aint(j));
        
        % Calcula el �rea de transferencia de masa en la partici�n actual de la torre
        % Primero halla el n�mero de gotas que hay en la partici�n
        N_g(j) = FV_E(j)*VP/V_g; 
        % Ahora calcula el �rea de TdeM en la partici�n
        A_M(j) = N_g(j)*As_g;
        
        % Resuelve el balance de masa por componente en Estado Estacionario, para 
        % hallar la conc. de Nicotina en la corriente de salida de Refinado de esta partici�n.
        X_A(j) = (X_A_EntPart(j)*dotn_RS - A_M(j)*N_A(j))/dotn_RS;
        
        % Resuelve el balance de masa por componente en EE para la concentraci�n de
        % Nicotina en la corriente de salidaentrada a la partici�n actual, que es igual a la
        % corriente de salida de la partici�n i+1, debajo de la partici�n actual (ver dibujo),
        % pues se est� operando a contracorriente
        Y_A(j+1) = (Y_A(j)*dotn_ES - A_M(j)*N_A(j))/dotn_ES;
        
        % Calcula el error de estimaci�n para decidir sobre la iteraci�n
        errorEst = abs(X_A(j)-X_Asup(j));
            
        % Ajusta el contador para tener registro del n�mero de iteraciones
        contador=contador+1;
        % Hace que el valor hallado en esta iteraci�n sea semilla por si se requiere otra iteraci�n
        X_Asup(j) = X_A(j);
    end
end

% Adiciona la concentraci�n de la entrada en el Refinado al vector de concentraciones
% hallado para tener toda la l�nea de operaci�n para graficarla y guardarla como el EE
X_A_EE = [X_A_Ent X_A];
% El vector de la Y ya est� correcto. Lo guarda con su nombre de EE
Y_A_EE = Y_A;
% El vector de �reas de transferencia calculadas en cada partici�n para este EE.
A_M_EE = A_M;
% El vector de fluxes de transferencias de masa calculadas
N_A_EE = N_A;

% Salva en un archivo el estado estacionario para hacer las simulaciones din�micas
save('EE_Torre.mat','L_T','D_T','dL','X_A_EE','Y_A_EE','A_M_EE','N_A_EE');

% Genera la curva de distribuci�n de equilibrio CDE para graficarla junto con la Ldeo
for l=1:600
    X_A_CDE(l) = (l-1)*0.00002;
    Y_A_CDE(l) = 16.07781726*X_A_CDE(l) +...
                (0.08593565238*X_A_CDE(l))/(X_A_CDE(l) + 0.0003360234205);
end

% Genera los puntos del equilibrio, que son los interceptos de las diferentes l�neas de
% fuerza impulsor con la CDE. Todo para graficarlos.
for j=1:NP
    Y_Aint(j) = 16.07781726*X_Aint(j) +...
                (0.08593565238*X_Aint(j))/(X_Aint(j) + 0.0003360234205);
end


figure(1), plot(X_A_EE,Y_A_EE,'or'), grid on, hold on, plot(X_LdeO,Y_LdeO,'.b'), 
            title('LdeO de Dise�o: ptos Azules y Real: circ. Rojos. CDE: circ.s negros. Interceptos: ptos verdes'),
            xlabel('X: Conc. en el Refinado [kmol Nicotina/kmol Agua]'), 
            ylabel('Y: Conc. en el Extracto [kmol Nicotina/kmol Tolueno]'),
            plot(X_A_CDE,Y_A_CDE,'ok'), plot (X_Aint,Y_Aint,'.g')
