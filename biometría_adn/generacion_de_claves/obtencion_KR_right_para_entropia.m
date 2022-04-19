function KR_right = obtencion_KR_right_para_entropia(V,paso,lambda,delta,aR,CR,KB1,T)

        % con este for se obtiene el punto mas alto de la funcion
        j = 1;
        for x = -150:paso:150
            f_R_right(j) = CR*((1./(1+exp(x))).*(exp(-((x-((lambda-delta+aR*V)/(KB1*T))).^2)*((KB1*T)/(4*lambda)))));
            
            if j == 1 
                punto_mas_alto = f_R_right(j);
            else
                if f_R_right(j) >= f_R_right(j-1)
                    punto_mas_alto = f_R_right(j);
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
            f_R_right(j) = CR*((1./(1+exp(x))).*(exp(-((x-((lambda-delta+aR*V)/(KB1*T))).^2)*((KB1*T)/(4*lambda)))));
            
            if f_R_right(j) >= punto_mas_alto_1_porct
                integral = integral + f_R_right(j);
            else
                integral = integral + 0;
            end
            
            j = j + 1;
        end              
        
        KR_right = integral*paso; 

end

