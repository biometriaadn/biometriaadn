clear all;
close all;
clc;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETROS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

paso = 0.01;
resolucion = 0.01;
limite_func = 1;
muestras = 20000000;

%% CALCULO ENTROPIA
% H = SUMA ( pi*Log2(pi))
i = 1;
for j = 0:1:muestras
%calculo de distribucion uniforme
% a = 9;
% b = 20;
% r = (b-a).*rand(1000,1) + a;

% N -- 9 a 20
% lambda -- 0.5 a 0.6
% aL -- 0.009 a 0.135
% aR -- 0.138 a 0.2 
% CL -- 1.7 a 7
% CR -- 0.7 a 3.5

Ni = 3;
Nf = 20;
N_al_1(i) = (Nf-Ni)*rand() + Ni;
N_al(i) = round(N_al_1(i));

lambdai = 0.5;
lambdaf = 0.6;
lambda_al(i) = (lambdaf-lambdai)*rand() + lambdai;

aLi = 0.009;
aLf = 0.135;
aL_al(i) = (aLf-aLi)*rand() + aLi;

aRi = 0.138;
aRf = 0.2;
aR_al(i) = (aRf-aRi)*rand() + aRi;

CLi = 1.7;
CLf = 7;
CL_al(i) = (CLf-CLi)*rand() + CLi;

CRi = 0.7;
CRf = 3.5;
CR_al(i) = (CRf-CRi)*rand() + CRi;

Vi = 2;
Vf = 10;
V_al(i) = (Vf-Vi)*rand() + Vi;

Ti = 273.15 - 50;
Tf = 273.15 + 50;
T_al(i) = (Tf-Ti)*rand() + Ti;

ki = 10^7;
kf = 10^9;
k_al(i) = (kf-ki)*rand() + ki;

var(i) = N_al(i)*lambda_al(i)*aL_al(i)*aR_al(i)*CL_al(i)*CR_al(i)*V_al(i);

i = i + 1;
end

i = 1;
for j = 0:1:muestras
    
    I_fun_ent(i) = obtencion_fun_ent(resolucion,paso,N_al(i),lambda_al(i),aL_al(i),aR_al(i),CR_al(i),CL_al(i),T_al(i),V_al(i),k_al(i));
    I_fun_ent_log(i) = log((I_fun_ent(i)));
    i = i + 1;
    
end

I_fun_ent_log_ordenado = sort(I_fun_ent_log);


for j = 1:1:muestras+1 
  
elementos(j) = 0; 
    
end

i = 1;
for j = 1:1:muestras+1
    if j == 1
        elementos(i) = elementos(i) + 1;
    else
    	if I_fun_ent_log_ordenado(j) == I_fun_ent_log_ordenado(j-1)
            elementos(i) = elementos(i) + 1;
        else
            i = i + 1;
            elementos(i) = elementos(i) + 1;
        end
    end
end

suma_prob = 0;
for j = 1:1:muestras+1
    
    probabilidad(j) = elementos(i)/(muestras+1);
    suma_prob = suma_prob + probabilidad(j);
    
end

H = 0;
for j = 1:1:muestras+1
    
    H = H + probabilidad(j)*log2(1/probabilidad(j));
    
end

suma_prob
H


figure
histogram(I_fun_ent)

figure
h = histogram(I_fun_ent_log,'Normalization','probability');
