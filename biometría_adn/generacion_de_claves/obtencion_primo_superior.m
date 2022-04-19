clear all;
close all;
clc;

%% obtencion numeros primos mas cercanos

N1=1659864707523991;
N2=1544077417120095;

es_primo = 0;
while (es_primo ~= 1)
    
    es_primo = isprime(N1);
    if (es_primo ~= 1)
        N1 = N1 + 1;
    end
    
end

es_primo = 0;
while (es_primo ~= 1)
    
    es_primo = isprime(N2);
    if (es_primo ~= 1)
        N2 = N2 + 1;
    end
    
end

numero_primo1=N1;
numero_primo2=N2;


















