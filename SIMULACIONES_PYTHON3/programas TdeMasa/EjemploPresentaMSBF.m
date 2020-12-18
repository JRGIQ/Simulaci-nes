%Ejemplo de la presentación de los MSBF notas de masa

close all
clear all
clc


%%%%%%%%%% Parámetros de operación del proceso
% Área superficial expuesta por las gotas
AM=140; %m2
% Masa cargada de gas (Refinado) en el equipo
MR(1)=5;     %kg de gas
% Masa cargada de Líquido (Extracto) en el equipo
ME(1)=20;   %kg de líquido
% Concentraciones iniciales de las fases en fracciones másicas libres Sto
Yw(1)=0.05;   %Fracción másica gas (Refinado)
Xw(1)=0.0095; %Fracción másica Líquido (Extracto)

% Con los datos disponibles se calcula la cantidad de solvente puro que existe
% en cada fase
% Primero se cambian las unidades de la concentración de los solventes de fracción
% másica respecto de solvente puro a fracción masa respecto de la fase como un todo 
w_Rs_R=1/(1+Yw(1));
w_Es_E=1/(1+Xw(1));
% Se calcula la cantidad de solvente que se cargo en cada fase: Rs y Es
% Solvente cargado con la fase gas (Refinado)
MRs=w_Rs_R*MR(1);
% Solvente cargado con la fase Líquida (Extracto)
MEs=w_Es_E*ME(1);

% Constantes del modelo
% Masa molecular del Soluto
MMA=14; %kg/kmol

%%%%%%%%%%% Datos para la simulación
TiempoSimula=20*60;   %minutos convertidos a segundos
Paso=1;  %Segundos
NumCiclos=ceil(TiempoSimula/Paso);


%%%%%%%%%% Ciclo for de solución del modelo usando Euler
for i=1:NumCiclos
    % Calcula la TdeM
    % Dato de concentración ficticia leída de la CDE
    XwAsterisco(i)=0.036-0.25*Xw(i);
    %Coeficiente englobante de TdeM para usar con conc. de líquido. Se
    %asume constante pero puede ser calculado con los kx y ky locales por
    %fase.
    KXw=1.3e-5;
    %Flux de TdeM usando conc. ficticia y valores del líquido 
    NA=KXw*(XwAsterisco(i)-Xw(i));
    
    % Solución de las Ecuaciones Diferenciales del modelo usando el método de Euler
    % Ecuación diferencial de la masa de Refinado en el sistema
    deltaMR=-AM*NA*MMA;
    MR(i+1)=MR(i)+deltaMR*Paso;

    % Ecuación diferencial de la concentración en el gas
    deltaYw=-(1/MRs)*(AM*NA*MMA);
    Yw(i+1)=Yw(i)+deltaYw*Paso;
    % Condicional para evitar concentraciones negativas en el Refinado
    if Yw(i+1)<0
        Yw(i+1)=0;
    end
    
    % Ecuación diferencial de la masa de Extracto en el sistema
    deltaME=AM*NA*MMA;
    ME(i+1)=ME(i)+deltaME*Paso;

    % Ecuación diferencial de la concentración en el gas
    deltaXw=(1/MRs)*(AM*NA*MMA);
    Xw(i+1)=Xw(i)+deltaXw*Paso;
    
end

figure(1), plot(MR,'*b'), grid on, title('Masa total del Refinado (gas)')
pause

figure(2), plot(Yw,'ob'), grid on, title('Concentración en frac. libre Solutode A en Refinado (gas)')
pause

figure(3), plot(ME,'*g'), grid on, title('Masa total del EXTRACTO (Liquido)')
pause

figure(4), plot(Xw,'og'), grid on, title('Conc. en frac. libre Sto de A en Extracto (Líquido) VERDE y Xw* NEGRO'), hold on,
plot (XwAsterisco,'ok')
