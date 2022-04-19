function I_fun_ent = obtencion_fun_ent(resolucion,paso,N_al,lambda_al,aL_al,aR_al,CR_al,CL_al,T_al,V_al,k_al)
%asociados a los 30nm
N = N_al;                        %adimensional
lambda = lambda_al;              %eV
delta = 0.55;                    %eV
aL = aL_al;                      %adimensional
aR = aR_al;                      %adimensional
aM = 1 - aL - aR;                %adimensional
CL = CL_al;                      %10^8 s^-1
CR = CR_al;                      %10^8 s^-1
V = V_al;

%generales
e = -1.6 * 10^-19;             %coulombs
KB = 1.380649 * (10^-23);      %(boltzmann constant) eV * K^-1                                  
KB1 = 8.6173324 *10^-5;        %(boltzmann constant) eV * K^-1
T = T_al;                 %grados kelvin (0º centigrados)
k = k_al;                      %10^8 s^-1



%% MODULO 2 --> Obtenemos Kf, Kb e integramos las funciones para obtener 
%               K_L_right, K_L_left, K_R_right, K_R_left

% obtenemos y representamos Kf 
Kf = k*(exp((aM*e*V)/(2*(N-1)*KB*T)));

% obtenemos y representamos Kb  
Kb = k*(exp((-aM*e*V)/(2*(N-1)*KB*T)));

% obtenemos y representamos KL_right
KL_right = obtencion_KL_right_para_entropia(resolucion,paso,lambda,delta,aL,CL,KB1,T);

% obtenemos y representamos KL_left
KL_left = obtencion_KL_left_para_entropia(resolucion,paso,lambda,delta,aL,CL,KB1,T);

% obtenemos y representamos KR_right
KR_right = obtencion_KR_right_para_entropia(resolucion,paso,lambda,delta,aR,CR,KB1,T);

% obtenemos y representamos KR_left
KR_left = obtencion_KR_left_para_entropia(resolucion,paso,lambda,delta,aR,CR,KB1,T);

%% MODULO 3 --> Obtenemos sistema de ecuaciones de puntos incoherentes
%               y la intensidad, y la representamos

%obtenemos puntos incoherentes
AA = KL_left + Kf;
BB = Kf * Kb;
CC = Kb + Kf;
DD = Kf * KL_right;

%inicializamos variables
for i = 1:1:N
    A(i) = 1;
    B(i) = 1;
    C(i) = 1;
    D(i) = 1;
end

%DEN
for i = 1:1:N
    if i == 1
        C(i) = AA;
    else
        C(i) = CC-(BB/(C(i-1)));
    end
end

%NUM_DEN
for i = 1:1:N
    if i == 1
        B(i) = 1;
    else
        B(i) = C(i-1);
    end
end

%NUM_NUM
for i = 1:1:N
    if i == 1
        A(i) = KL_right;
    elseif i == 2
        A(i) = DD;
    else
        A(i) = A(i-1)/C(i-2);
    end
end

%FUNCION
for i = N:-1:1
    if i == N
        D_n(i) = ((A(i)/B(i))+(Kb))/(C(i));
    else
        D_n(i) = ((A(i)/B(i))+(Kb*D(i+1)))/(C(i));
    end
end

%FUN PLR
for i = 1:1:N
    PLR_var = 1+D(i); 
end
PLR = 1/PLR_var;

for i = N:-1:1
    if i == N
        P(i) = (((A(i)/B(i))*PLR)+(Kb*PLR))/(C(i));
    else
        P(i) = (((A(i)/B(i))*PLR)+(Kb*P(i+1)))/(C(i));
    end
end

%obtenemos I
I_fun_ent = e * (KL_right*PLR - KL_left*P(1));
    

end