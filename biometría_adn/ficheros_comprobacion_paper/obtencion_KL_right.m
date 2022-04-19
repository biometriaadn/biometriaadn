function KL_right = obtencion_KL_right(resolucion,paso,lambda,delta,aL,CL,KB1,T,V_min,V_max)

    % recorremos con un for el voltaje 
    i = 1;
    for V = V_min:resolucion:V_max
        
        % con este for se obtiene el punto mas alto de la funcion
        j = 1;
        for x = -150:paso:150
            f_L_right(j) = CL*((1./(1+exp(x))).*(exp(-((x-((lambda+delta+aL*V)/(KB1*T))).^2)*((KB1*T)/(4*lambda)))));
            
            if j == 1 
                punto_mas_alto = f_L_right(j);
            else
                if f_L_right(j) >= f_L_right(j-1)
                    punto_mas_alto = f_L_right(j);
                end
            end
            
            j = j + 1;
            
        end        
        
        % multiplicamos por 0.01 el punto mas alto para obtener los limites
        punto_mas_alto_1_porct = punto_mas_alto*0.01;
        
        
        % con este for se obtiene la integral dela funcion y si el valor 
        % de la funcion es menor del 1% dl valor del dato mas alto lo 
        % tomamos como 0
        j = 1;
        integral = 0;
        for x = -150:paso:150
            f_L_right(j) = CL*((1./(1+exp(x))).*(exp(-((x-((lambda+delta+aL*V)/(KB1*T))).^2)*((KB1*T)/(4*lambda)))));
            
            if f_L_right(j) >= punto_mas_alto_1_porct
                integral = integral + f_L_right(j);
            else
                integral = integral + 0;
            end
            
            j = j + 1;
        end              
        
        KL_right(i) = integral*paso; 
        
        i = i + 1;
        
    end 
end