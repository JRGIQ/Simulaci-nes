%%%%%%%% IDENTIFICACACION DEL PROCESO
% Este programa permite hallar el estado estacionario de un equipo de
% deshumidificacion
% operando en Contra-Corriente. Dicho estado estacionario se puede usar solo en este proceso.

% Autores: JULIANA- DAVID-JOSE
% Última actualización: Abri, 2018

% Nomenclatura usada y esquemático de corrientes:
% R: Refinado. SdeP I. Solucion de aire con benceno. R solución
% completa, RS solvente en la solución de la fase R
% E: Extracto. SdeP II. Benceno liquido(GOTAS). E Benceno liquido
% X_A: Concentracion de benceno en el Extracto [kmol benceno/kmol
% benceno]que va hacer igual a uno  por tener benceno como soluto y
% reefinado
% Y_A: Concentración de A la fase Refinado [kmol benceno /kmol aire]

%    X_A(0)      Y_A(1)    
%       -----------    ---- Plano 1, partición 1:  Entrada del Extracto (Salida de R)
%       |         |
%       |         |
%       |  Torre  |
%       |         |
%       |         |
%       -----------    ---- Plano 2, partición NP:  Entrada del Refinado (Salida de E)
%    X_A(NP)   Y_A(NP+1) 

% En términos de vectores, la posición j siempre indica flujos y concentraciones correspondientes
% a la partición j, pero en sus salidas, NO en las entradas. Por eso, al graficar X(1),Y(1), 
% no se verá la entrada del Refinado a la torre sino la salida de R de la primera partición
% con la salida de la primera partición o salida de E de la torre. Debe prestarse atención
% a esto porque estas parejas forman la Línea de Operación, pero no son las que se usan para el
% cálculo de las propiedades dentro de la partición, pues una es salida y la otra es entrada.

%%%%%%%%% AJUSTE DE LA PLATAFORMA DE PROGRAMACIÓN
clear all 
close all
clc
format long

%%%%%%%% INICIALIZACIÓN DE VARIABLES GENERALES
% Altura de la partición a tomar en la torre
dL = 0.05;                          % [m] 
% Altura útil de la torre
L_T = 3.1;                   %[m]
% Número de particiones hechas en la torre para resolver el modelo
NP = ceil(L_T/dL);       % [Adim] 
% Diámetro de la torre
D_T = 1;                        %[m]
% Tolerancia del error en la iteración para convergencia
Tol_error = 1e-8;%
% Valor máximo de ciclos de iteración de convergencia
MaxContador = 100;

%%%%%%%% DECLARACIÓN DE CONSTANTES DEL MODELO
% Aceleración de la gravedad
g = 9.8;                        %[m/s2]
% Tensión superficial de la fase dispersa: solución de Tolueno
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
% de operar a otra T se deberán ajustar a los valores correctos.
% Densidad del Tolueno a 303K
rho_T = 857.7;                   %[kg/m3]
% Densidad del Agua a 303K
rho_W = 995.7;                   %[kg/m3]
% Densidad de la Nicotina a 303K
rho_N = 1010;                    %[kg/m3]

%%%%%%%%% INICIALIZACIÓN DE VARIABLES ESCALARES DEL MODELO
% Diámetro de la gota de la fase dispersa: solución de Nicotina en Tolueno,
% dada en milímetros pero convertida a metros.
D_g = 1/1000;                      % [m]
% Flujo molar de solvente en el Refinado: Agua, que es cte en toda la torre
% dado en moles por segundo pero convertido a kmoles
dotn_RS = 0.2;                % [kmol/s] de Agua
% Flujo molar de se solvente en el Extracto: Tolueno, cte en toda la torre
% dado en moles por segundo pero convertido a kmoles
dotn_ES = 0.0158;            % [kmol/s] de Tolueno
% Concentraciones en el Plano 1, Superior:
% Concentración de Entrada (dato) del Refinado. Corresponde a la posición Cero 0
X_A_Ent = 0.0087;                %[kmol Nicotina/kmol RS]
% Concentración de Salida Esperada del Extracto
Y_A(1) = 0.1;                 %[kmol Nicotina/kmol ES]
% Concentraciones en el Plano 2, Inferior:
% Concentración de Entrada (dato) del Extracto. Corresponde a la posición NP+1
Y_A_Ent = 0.005;              %[kmol Nicotina/kmol ES]
% Concentración de Salida Esperada del Refinado
X_A(NP) = 0.0012;               %[kmol Nicotina/kmol RS]

% Cálculos usando las inicializaciones previas.
% Área de flujo total a torre vacía
A_FT = pi*(D_T^2)/4;            %[m2]
% Volumen de una partición
VP = A_FT*dL;                   %[m3]
% Intervalo total de concentraciones del Refinado (con la Conc. de salida Esperada)
Int_X_A = X_A_Ent - X_A(NP);    %[kmol Nicotina/kmol RS]
% Calcula la pendiente de la línea de operación que se espera de la torre
m_LdeO = (Y_A_Ent - Y_A(1))/(X_A(NP) - X_A_Ent);    %[kmol RS/kmol ES]
% Volumen de una gota
V_g = (1/6)*pi*D_g^3;           %[m3]
% Área superficial de una gota
As_g = pi*D_g^2;                %[m2]

%%%%%%%%% INICIALIZACIÓN DE VARIABLES VECTORIALES DEL MODELO
% Genera puntos para graficar la línea de operación LdeO. Asigna la primera
% pareja, correspondiente a las concentraciones en el Plano 1
X_LdeO(1)= X_A_Ent;
Y_LdeO(1) = Y_A(1);
% Ciclo for para generar el resto de puntos de la línea de operación LdeO
for l = 1:NP
    % Valor de concentración supuesta para el Refinado: con espaciado uniforme, a la entrada
    % de la partición l+1 que es lo mismo que la salida de la partición l.
    X_LdeO(l+1) = X_LdeO(l)-(Int_X_A/NP);
    % Calcula el valor correspondiente de la concentración de Extracto, que es la de salida de
    % la partición l+1, usando la ecuación de la línea de operación
    Y_LdeO(l+1) = m_LdeO*(X_LdeO(l+1) - X_A(NP)) + Y_A_Ent;
    % Nótese que la LdeO contiene siempre la coordenada X de la entrada a la partición, 
    % mientras que la coordenada Y corresponde a la salida de esa partición. Es decir,
    % no son puntos de salidas de la partición tratada.
end

% Inicializa el vector de concentraciones finales a obtener del Refinado, en el cual
% la primera posición es conocida e igual a la concentración de entrada de R dada
X_A(1) = X_A_Ent;

%%%%%%%%%% RESUELVE TODA LA TORRE POR PARTICIONES HASTA LOGRAR LA CONVERGENCIA
for j = 1:NP  % Ciclo for para la variación de la Posición de la Partición en la Torre
    % Asigna una concentración semilla para la iteración en la salida del Refinado
    % de la partición actual, usando una repartición uniforme de distribución de 
    % concentración. Parte de la concentración conocida de la entrada de R.
    X_Asup(j) = X_A_Ent - j*(Int_X_A/NP);
    
    if j==1
        % Asigna la concentración de Refinado en la corriente de entrada a la partición
        % para evitar el problema del índice cero en Matlab para los vectores
        X_A_EntPart(j) = X_A_Ent; % La primera vez es la de entrada
    else
        % Asigna la concentración de Refinado en la corriente de entrada a la partición
        % para el resto de particiones, que no presentan el problema del índice cero en Matlab
        X_A_EntPart(j) = X_A(j-1); % Las otras veces es la anterior
    end    
    
    % La concentración de la fase Extracto a la salida de esta partición se conoce para
    % la primera partición, pues es la concentración de salida esperada. Para las demás
    % particiones será un valor que la iteración previa ya calculó como el valor 
    % definitivo del estado estacionario.
    if j==1
        Y_A(j) = Y_A(1);
    end
          
    % Ciclo while de iteración para buscar convergencia en la concentración del 
    % Refinado a la salida de cada partición del equipo.
    % Inicializa el error en un valor dos veces la torelacia para entrar a la iteración
    errorEst = 2*Tol_error;
    % Inicializa contador para el número de iteraciones máximas
    contador = 0; 
    while errorEst > Tol_error && contador < MaxContador
        % Bloque de cálculo del Flux de trasferencia de masa.
        % Como se requiere tener las concentraciones en las unidades adecuadas para
        % actualizar las propiedadess de las fases, hace las conversiones necesarias.
        % Conversión a fracción molar respecto de la fase como un todo
        x_Asup(j) = X_Asup(j)/(1+X_Asup(j));
        y_A(j) = Y_A(j)/(1+Y_A(j)); 
        
        % Conversión a fracción másica respecto del solvente puro
        Xw_Asup(j) = X_Asup(j)*MM_N/MM_W;
        Yw_A(j) = Y_A(j)*MM_N/MM_T;
    
        % Conversión a fracción másica respecto de la fase como un todo
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
        
        % Como va a requerir flujos volumétricos de las fases para los
        % cálculos, los halla con los datos de flujos molares existentes.
        % Calcula el flujo total de Refinado como la suma del solvente y el soluto 
        dotn_R(j) = dotn_RS + dotn_RS*X_Asup(j);        %[kmol/s]
        % Lo convierte a flujo másico con la masa molecular del Refinado en ese punto
        dotm_R(j) = dotn_R(j)*MM_R(j);                  %[kg/s]
        % Ahora lo conviere a flujo volumétrico usando la densidad actual del Refinado
        dotV_R(j) = dotm_R(j)/rho_R(j);                 %[m3/s]
        
        % Calcula el flujo total de Extracto como la suma del solvente y el soluto 
        dotn_E(j) = dotn_ES + dotn_ES*Y_A(j);           %[kmol/s]
        % Lo convierte a flujo másico con la masa molecular del Extracto en ese punto
        dotm_E(j) = dotn_E(j)*MM_E(j);                  %[kg/s]
        % Ahora lo conviere a flujo volumétrico usando la densidad actual del Refinado
        dotV_E(j) = dotm_E(j)/rho_E(j);                 %[m3/s]
        
        % Calcula el coeficiente local de TdeM para la Fase Refinado, con concentraciones X
        % Cálculo del número de Reynolds de la gota, pero primero halla todo lo previo
        % Fracción volumétrica del Extracto o hold-up, fase dispersa, en la torre
        FV_E(j) = dotV_E(j)/(dotV_R(j) + dotV_E(j));
        % Velocidad de la fase Refinado en la torre
        v_R(j) = dotV_R(j)/(A_FT*(1-FV_E(j)));          %[m/s]
        % Velocidad de la fase Extracto en la torre
        v_E(j) = dotV_E(j)/(A_FT*FV_E(j));              %[m/s]
        % Velocidad de deslizamiento entre las fases.
        v_s(j) = v_R(j) + v_E(j);                       %[m/s]
        % Número de Reynolds de la gota, porque así lo pide la correlación
        Re_g(j) = ((D_g*v_s(j)*rho_R(j))/mu_R(j));
        
        % Cálculo del número de Schmidt
        Sc_R(j) = mu_R(j)/(rho_R(j)*Dif_N_W);  
        
        % Coeficiente propio de la correlación de Skellan
        K_Skll = (Re_g(j)^(1/8))*((mu_R(j)/mu_E(j))^(1/4))*((mu_R(j)*v_s(j)/sigma_T)^(1/6));
        % Factor de corrección por forma
        F_Skll = 0.281+ 1.615*K_Skll+3.73*K_Skll^2-1.874*K_Skll^3;
        % Coeficiente local de TdeM para la fase continua, el Refinado
        k_X(j) = (Dif_N_T/D_g)*(2+0.463*Re_g(j)^(0.484)*Sc_R(j)^(0.339)*...
                    ((D_g*g^(1/3))/(Dif_N_W^(2/3)))^(0.072))*F_Skll;
        
        % Calcula coeficiente local de TdeM  para Fase Extracto, con concentraciones Y
        k_Y(j) = 0.00375*v_s(j)/(1+(mu_E(j)/mu_R(j)));
        
        % Calcula la concentración en la interfase, para lo cual primero halla la
        % pendiente de la línea de Fuerza Impulsora (FI)
        m_FI = -k_Y(j)/k_X(j);
        % Resuelve la ecuación resultante del cruce de la FI y la CDE usando el punto actual
        % Usa la función que encuentra el cero o la raíz de una función no lineal
        % Declara la función Cruce de FI y CDE en un .m externo y declara una función
        % muda FunConArg que llama a la función original para pasarle los parámetros
        FunConArg=@(x)CruceFI_CDE(x,X_Asup(j),Y_A(j),m_FI);
        % Da un punto inicial para buscar la solución
        if j==1
            % Para la primera partición da un valor arbitrario que esté cerca del
            % valor esperado. Este valor se determinan mirando la LdeO y la CDE.
            x_0 = 0.001;
        else
            % Para las demás particiones toma como semilla el valor anterior
            % buscando mantener la continuidad de la salida 
            x_0 = X_Aint(j-1);
        end
        % Llama la función que encuntra la raíz o el cero de la función no lineal
        X_Aint(j) = fzero(FunConArg,x_0);
        
        % Calcula el flux de TdeM usando la ecuación con coeficiente local y las
        % concentraciones asociadas con la fase R
        N_A(j) = k_X(j)*(X_Asup(j) - X_Aint(j));
        
        % Calcula el área de transferencia de masa en la partición actual de la torre
        % Primero halla el número de gotas que hay en la partición
        N_g(j) = FV_E(j)*VP/V_g; 
        % Ahora calcula el área de TdeM en la partición
        A_M(j) = N_g(j)*As_g;
        
        % Resuelve el balance de masa por componente en Estado Estacionario, para 
        % hallar la conc. de Nicotina en la corriente de salida de Refinado de esta partición.
        X_A(j) = (X_A_EntPart(j)*dotn_RS - A_M(j)*N_A(j))/dotn_RS;
        
        % Resuelve el balance de masa por componente en EE para la concentración de
        % Nicotina en la corriente de salidaentrada a la partición actual, que es igual a la
        % corriente de salida de la partición i+1, debajo de la partición actual (ver dibujo),
        % pues se está operando a contracorriente
        Y_A(j+1) = (Y_A(j)*dotn_ES - A_M(j)*N_A(j))/dotn_ES;
        
        % Calcula el error de estimación para decidir sobre la iteración
        errorEst = abs(X_A(j)-X_Asup(j));
            
        % Ajusta el contador para tener registro del número de iteraciones
        contador=contador+1;
        % Hace que el valor hallado en esta iteración sea semilla por si se requiere otra iteración
        X_Asup(j) = X_A(j);
    end
end

% Adiciona la concentración de la entrada en el Refinado al vector de concentraciones
% hallado para tener toda la línea de operación para graficarla y guardarla como el EE
X_A_EE = [X_A_Ent X_A];
% El vector de la Y ya está correcto. Lo guarda con su nombre de EE
Y_A_EE = Y_A;
% El vector de áreas de transferencia calculadas en cada partición para este EE.
A_M_EE = A_M;
% El vector de fluxes de transferencias de masa calculadas
N_A_EE = N_A;

% Salva en un archivo el estado estacionario para hacer las simulaciones dinámicas
save('EE_Torre.mat','L_T','D_T','dL','X_A_EE','Y_A_EE','A_M_EE','N_A_EE');

% Genera la curva de distribución de equilibrio CDE para graficarla junto con la Ldeo
for l=1:600
    X_A_CDE(l) = (l-1)*0.00002;
    Y_A_CDE(l) = 16.07781726*X_A_CDE(l) +...
                (0.08593565238*X_A_CDE(l))/(X_A_CDE(l) + 0.0003360234205);
end

% Genera los puntos del equilibrio, que son los interceptos de las diferentes líneas de
% fuerza impulsor con la CDE. Todo para graficarlos.
for j=1:NP
    Y_Aint(j) = 16.07781726*X_Aint(j) +...
                (0.08593565238*X_Aint(j))/(X_Aint(j) + 0.0003360234205);
end


figure(1), plot(X_A_EE,Y_A_EE,'or'), grid on, hold on, plot(X_LdeO,Y_LdeO,'.b'), 
            title('LdeO de Diseño: ptos Azules y Real: circ. Rojos. CDE: circ.s negros. Interceptos: ptos verdes'),
            xlabel('X: Conc. en el Refinado [kmol Nicotina/kmol Agua]'), 
            ylabel('Y: Conc. en el Extracto [kmol Nicotina/kmol Tolueno]'),
            plot(X_A_CDE,Y_A_CDE,'ok'), plot (X_Aint,Y_Aint,'.g')
