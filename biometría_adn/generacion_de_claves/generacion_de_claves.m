clear all;
close all;
clc;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETROS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

paso = 0.1;
resolucion = 0.1;
num_v = 2;
limite_func = 1;
muestras = 1000000;
redondeo_el = 15;

%% GENERACION CLAVES

%determinamos los parametros de forma aleatoria
i = 1;
for j = 1:1:muestras
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

for t = 1:1:num_v
    Vi = 2;
    Vf = 10;
    V_al(t,i) = (Vf-Vi)*rand() + Vi;
end

Ti = 273.15 - 50;
Tf = 273.15 + 50;
T_al(i) = (Tf-Ti)*rand() + Ti;

ki = 10^7;
kf = 10^9;
k_al(i) = (kf-ki)*rand() + ki;

i = i + 1;
end

%obtenemos la clave, y su logaritmo(que le hace aumentar la entropia)
for j = 1:1:muestras
 for i = 1:1:num_v
    I_fun_clave(i,j) = obtencion_fun_ent(resolucion,paso,N_al(j),lambda_al(j),aL_al(j),aR_al(j),CR_al(j),CL_al(j),T_al(j),V_al(i,j),k_al(j));
    I_fun_clave_log(i,j) = log((I_fun_clave(i,j)));
 end
end

%redondeamos y quitamos los decimales a las claves
for j = 1:1:muestras
    for i = 1:1:num_v
        I_fun_clave_redondeada(i,j) = round(abs(I_fun_clave_log(i,j)*10^redondeo_el));
    end
end


%pasamos las claves a binario y juntamos los "cachos"
filas = 1;
columnas = 1;
for j = 1:1:muestras
    var_cambio = 0;
    for i = 1:1:num_v
        var_cambio = var_cambio + 1;
        if I_fun_clave_log(i,j) < 0 
            claves_binaria(filas,columnas) = int8(1);
            columnas = columnas + 1;
        else
            claves_binaria(filas,columnas) = int8(0);
            columnas = columnas + 1;
        end

        I_fun_clave_redondeada_abs = abs(I_fun_clave_redondeada(i,j));
        %codigo que transforma elnumero en binario una vez redondeado, etc
        %desde aqui--------------------------------------------------------
        var_num_dec = I_fun_clave_redondeada_abs;
        var_num_div = var_num_dec;
        var_num_res = var_num_div;
        var_i = 1;
        while var_num_div ~= 1
            var_resto(var_i) = mod(var_num_div,2);
            var_num_div = fix(var_num_div/2);
            var_i = var_i + 1;

            if var_num_div == 1
                var_resto(var_i) = 1;
            end

        end

        var_down_j = length(var_resto);
        for var_j = 1:1:length(var_resto)

            num_binario(var_down_j) = var_resto(var_j);
            var_down_j = var_down_j - 1;

        end
        %hasta aqui--------------------------------------------------------
        for t = 1:1:length(num_binario)
            claves_binaria(filas,columnas) = int8(num_binario(1,t));
            columnas = columnas+1;
        end
        
        if var_cambio == num_v
            columnas = 1;
        end
    end
    filas = filas + 1;
end

%pasamos las claves a decimal
[filas_clave_binaria,columnas_clave_binaria] = size(claves_binaria);
var_down = columnas_clave_binaria;
var_t = 1;
claves_num(1)=0;
for var_t1 = 1:1:filas_clave_binaria
    num_temp = 0;
    for var_t2 = 1:1:columnas_clave_binaria
        if claves_binaria(var_t1,var_t2) == 1
            num_temp = num_temp + 2^var_down;
            var_down = var_down - 1;
        else
            var_down = var_down - 1; 
        end
    end
    claves_num(var_t) = num_temp;
    var_t = var_t + 1;
    var_down = columnas_clave_binaria;
end


%% Calculo de la entropia de las claves finales

claves_num_ordenado = sort(claves_num);

for j = 1:1:muestras 
  
elementos(j) = 0; 
    
end

i = 1;
for j = 1:1:muestras
    if j == 1
        elementos(i) = elementos(i) + 1;
    else
    	if claves_num_ordenado(j) == claves_num_ordenado(j-1)
            elementos(i) = elementos(i) + 1;
        else
            i = i + 1;
            elementos(i) = elementos(i) + 1;
        end
    end
end

suma_prob = 0;
for j = 1:1:muestras
    
    probabilidad(j) = elementos(i)/(muestras+1);
    suma_prob = suma_prob + probabilidad(j);
    
end


% H = SUMA ( pi*Log2(pi))
H = 0;
for j = 1:1:muestras
    
    H = H + probabilidad(j)*log2(1/probabilidad(j));
    
end

