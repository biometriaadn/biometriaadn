clear all;
close all;
clc;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETROS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

paso = 0.1;
resolucion = 0.1;
limite_func = 1;
muestras = 10;
V_min = -5;
V_max = 5;

%asociados a los 35nm
N = 10;                                         %adimensional
lambda = 0.6;                                   %eV
delta = 0.55;                                   %eV
aL = 0.135;                                     %adimensional
aR = 0.200;                                       %adimensional
aM = 1 - aL - aR;                               %adimensional
CL = 2 * 10^8;                                %10^8 s^-1
CR = 1 * 10^8;                                  %10^8 s^-1

%generales
e = -1.6 * 10^-19;                              %coulombs
KB = 1.380649 * (10^-23);                       %(boltzmann constant) eV * K^-1                                  
KB1 = 8.6173324 *10^-5;                         %(boltzmann constant) eV * K^-1
T = 273.15-50;                                  %grados kelvin (0º centigrados)
volt = [V_min:resolucion:V_max];                %usada para representar en funcion del voltaje
volt1 = [V_min:resolucion:V_max];               %usada para representar en funcion del voltaje
k = 1*10^8;                                       %10^8 s^-1
a= [-200:paso:200];                             %usada para representar

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODELO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MODULO 1 --> Representamos las funciones que se van a integrar
% 
% %funcion L_right
% j = 1;
% for V = V_min:resolucion:V_max
%     
%     for x = -150:paso:150
%         f_L_right(j) = CL*((1./(1+exp(x))).*(exp(-((x-((lambda+delta+aL*V)/(KB1*T))).^2)*((KB1*T)/(4*lambda)))));
%         
%         j = j + 1;
%     end
%     
%     figure
%     plot(ex,f_L_right)
%     ylabel('f_L_right') 
%     xlabel('x')
%     
%     j = 1;
%     
% end
% 
% %funcion L_left
% j = 1;
% for V = V_min:resolucion:V_max
%     
%     for x = -150:paso:150
%         f_L_left(j) = CL*((1./(1+exp(x))).*(exp(-((x-((lambda-delta-aL*V)/(KB1*T))).^2)*((KB1*T)/(4*lambda)))));
%         
%         j = j + 1;
%     end
%     
%     figure
%     plot(ex,f_L_left)
%     ylabel(['f_L_left: V = ',num2str(V)]) 
%     xlabel('x')
% 
%     j = 1;
%     
% end
% 
% %funcion R_right
% j = 1;
% for V = V_min:resolucion:V_max
%     
%     for x = -150:paso:150
%         f_R_right(j) = CR*((1./(1+exp(x))).*(exp(-((x-((lambda-delta+aR*V)/(KB1*T))).^2)*((KB1*T)/(4*lambda)))));
%     
%         j = j + 1;
%     end
%     
%     figure
%     plot(ex,f_R_right) 
%     ylabel(['f_R_right: V = ',num2str(V)])
%     xlabel('x')
% 
%     j = 1;
%     
% end
% 
% %funcion R_left
% j = 1;
% for V = V_min:resolucion:V_max
%     
%     for x = -150:paso:150
%         f_R_left(j) = CR*((1./(1+exp(x))).*(exp(-((x-((lambda+delta-aR*V)/(KB1*T))).^2)*((KB1*T)/(4*lambda)))));
%     
%         j = j + 1;
%     end
%     
%     figure
%     plot(ex,f_R_left)
%     ylabel(['f_R_left: V = ',num2str(V)]) 
%     xlabel('x')
%     
%     j = 1;
%     
% end


%% MODULO 2 --> Obtenemos Kf, Kb e integramos las funciones para obtener K_L_right, K_L_left, K_R_right, K_R_left.       

% obtenemos y representamos Kf 
j = 1;
for V = V_min:resolucion:V_max
    
    Kf(j) = k*(exp((aM*e*V)/(2*(N-1)*KB*T)));
    
    j = j + 1;
    
end

figure
subplot(2,3,3)
plot(volt,Kf)
ylabel('Kf') 
xlabel('V')

% obtenemos y representamos Kb 
j = 1;
for V = V_min:resolucion:V_max
    
    Kb(j) = k*(exp((-aM*e*V)/(2*(N-1)*KB*T)));
    j = j + 1;
    
end

% figure
subplot(2,3,6)
plot(volt,Kb)
ylabel('Kb') 
xlabel('V')

% obtenemos y representamos KL_right
KL_right = obtencion_KL_right(resolucion,paso,lambda,delta,aL,CL,KB1,T,V_min,V_max);

% figure
subplot(2,3,1)
plot(volt,KL_right)
ylabel('KL right') 
xlabel('V')

% obtenemos y representamos KL_left
KL_left = obtencion_KL_left(resolucion,paso,lambda,delta,aL,CL,KB1,T,V_min,V_max);

% figure
subplot(2,3,4)
plot(volt1,KL_left)
ylabel('KL left') 
xlabel('V')

% obtenemos y representamos KR_right
KR_right = obtencion_KR_right(resolucion,paso,lambda,delta,aR,CR,KB1,T,V_min,V_max);

% figure
subplot(2,3,2)
plot(volt,KR_right)
ylabel('KR right') 
xlabel('V')

% obtenemos y representamos KR_left
KR_left = obtencion_KR_left(resolucion,paso,lambda,delta,aR,CR,KB1,T,V_min,V_max);

% figure
subplot(2,3,5)
plot(volt,KR_left)
ylabel('KR left') 
xlabel('V')


%% MODULO 3 --> Obtenemos sistema de ecuaciones de puntos incoherentes y la intensidad, y la representamos              

j = 1;
for V = V_min:resolucion:V_max
    %obtenemos puntos incoherentes
    A = KL_left(j) + Kf(j);
    B = Kf(j) * Kb(j);
    C = Kb(j) + Kf(j);
    D = Kf(j) * KL_right(j);

    %DEN
    C10 = C-(B/C-(B/C-(B/(C-(B/(C-(B/(C-(B/(C-(B/(C-(B/(C-(B/A)))))))))))))));
    C9 = C-(B/C-(B/(C-(B/(C-(B/(C-(B/(C-(B/(C-(B/(C-(B/A))))))))))))));
    C8 = C-(B/(C-(B/(C-(B/(C-(B/(C-(B/(C-(B/(C-(B/A)))))))))))));
    C7 = C-(B/(C-(B/(C-(B/(C-(B/(C-(B/(C-(B/A)))))))))));
    C6 = C-(B/(C-(B/(C-(B/(C-(B/(C-(B/A)))))))));
    C5 = C-(B/(C-(B/(C-(B/(C-(B/A)))))));
    C4 = C-(B/(C-(B/(C-(B/A)))));
    C3 = C-(B/(C-(B/A)));
    C2 = C-(B/A);
    C1 = A;

    %NUM_DEN
    B10 = C9;
    B9 = C8;
    B8 = C7;
    B7 = C6;
    B6 = C5;
    B5 = C4;
    B4 = C3;
    B3 = C2;
    B2 = C1;
    B1 = 1;

    %NUM_NUM
    A10 = (((((((D/A)/C2)/C3)/C4)/C5)/C6)/C7)/C8;
    A9 = ((((((D/A)/C2)/C3)/C4)/C5)/C6)/C7;
    A8 = (((((D/A)/C2)/C3)/C4)/C5)/C6;
    A7 = ((((D/A)/C2)/C3)/C4)/C5;
    A6 = (((D/A)/C2)/C3)/C4;
    A5 = ((D/A)/C2)/C3;
    A4 = (D/A)/C2;
    A3 = D/A;
    A2 = D;
    A1 = KL_right(j);

    %FUNCION
    D10 = ((A10/B10)+(Kb(j)))/(C10);
    D9 = ((A9/B9)+(Kb(j)*D10))/(C9);
    D8 = ((A8/B8)+(Kb(j)*D9))/(C8);
    D7 = ((A7/B7)+(Kb(j)*D8))/(C7);
    D6 = ((A6/B6)+(Kb(j)*D7))/(C6);
    D5 = ((A5/B5)+(Kb(j)*D6))/(C5);
    D4 = ((A4/B4)+(Kb(j)*D5))/(C4);
    D3 = ((A3/B3)+(Kb(j)*D4))/(C3);
    D2 = ((A2/B2)+(Kb(j)*D3))/(C2);
    D1 = ((A1/B1)+(Kb(j)*D2))/(C1);

    %FUN PLR
    PLR(j) = 1/(1+D10+D9+D8+D7+D6+D5+D4+D3+D2+D1);
    P10(j) = (((A10/B10)*PLR(j))+(Kb(j)*PLR(j)))/(C10);
    P9(j) = (((A9/B9)*PLR(j))+(Kb(j)*PLR(j)))/(C9);
    P8(j) = (((A8/B8)*PLR(j))+(Kb(j)*P9(j)))/(C8);
    P7(j) = (((A7/B7)*PLR(j))+(Kb(j)*P8(j)))/(C7);
    P6(j) = (((A6/B6)*PLR(j))+(Kb(j)*P7(j)))/(C6);
    P5(j) = (((A5/B5)*PLR(j))+(Kb(j)*P6(j)))/(C5);
    P4(j) = (((A4/B4)*PLR(j))+(Kb(j)*P5(j)))/(C4);
    P3(j) = (((A3/B3)*PLR(j))+(Kb(j)*P4(j)))/(C3);
    P2(j) = (((A2/B2)*PLR(j))+(Kb(j)*P3(j)))/(C2);
    P1(j) = (((A1/B1)*PLR(j))+(Kb(j)*P2(j)))/(C1);

    %obtenemos I
    I1(j) = e * (KL_right(j)*PLR(j) - KL_left(j)*P1(j));
    j = j + 1;
end

%representamos los puntos incoherentes y la intensidad
figure
plot(volt,I1)
ylabel('I') 
xlabel('V')

% figure
% plot(volt,P1)
% ylabel('P1') 
% xlabel('V')
% 
% figure
% plot(volt,P2)
% ylabel('P2') 
% xlabel('V')
% 
% figure
% plot(volt,P3)
% ylabel('P3') 
% xlabel('V')
% 
% figure
% plot(volt,P4)
% ylabel('P4') 
% xlabel('V')
% 
% figure
% plot(volt,P5)
% ylabel('P5') 
% xlabel('V')
% 
% figure
% plot(volt,P6)
% ylabel('P6') 
% xlabel('V')
% 
% figure
% plot(volt,P7)
% ylabel('P7') 
% xlabel('V')
% 
% figure
% plot(volt,P8)
% ylabel('P8') 
% xlabel('V')
% 
% figure
% plot(volt,P9)
% ylabel('P9') 
% xlabel('V')
% 
% figure
% plot(volt,PLR)
% ylabel('PLR') 
% xlabel('V')
