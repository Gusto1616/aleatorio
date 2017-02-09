
clc; clear all; close all;

disp('Running...')

%Pathway Parameters:
%Distance and Time
xf=10000;%m
xi=0;%m
t4=600; %s
t0=0; %s
Travel_d=xf-xi; % distance in m
Travel_t=t4-t0; % total time in s

%Constantes para fazer a convers�o entre as unidades
kmh2ms=(1000/3600);
km2m=1000;
D_tolerance=5; %Duistance in meters
V_tolerance=1*kmh2ms;

state=0;

%Necess�rio primeiro ser definido a distancia e o tempo de percurso
lengthTime_vector=t4+1;
time_vector=ones(lengthTime_vector,1);

%Distance and Velocity vectors - Initialize with all 0
x_vector=zeros(lengthTime_vector,1);
v_vector=zeros(lengthTime_vector,1);

%Time Vector - Constru��o do tempo do vector discretizado de 1 em 1 segundo
time_vector=time_vector.*(t0:1:t4)';

%Chama Ficheiro .m com os parametros do comboio utilizado
TrainParameters

%Track Characteristics: Vmax and Slope
%Vmax=[0 50;1 120;4 100;7 90]; %Feito para receber as distancias em km e os limites em km/h
Vmax=[0 50;1 120;4 50;7 190];

%Starting with average value
Vavg=Travel_d/Travel_t;

%random value for acceleration
%a1=accGen(acc_min,acc_max);
a1=0.5;
a2=0;
a3=-0.055;
a4=-0.5;
Vop=90*kmh2ms;
Vi=0*kmh2ms; %velocidade inicial
Vf=0*kmh2ms;
Vavg=Travel_d/Travel_t;

limits_number=size(Vmax,1);
points_number=size(Vmax,1);
Vmax_vector=zeros(length(x_vector),1);

[t1_aux,t2_aux,t3_aux]=times_calculation(Vi,Vf,Vop,t0,t4,a1,a2,a3,a4,Travel_d);
Vbk=a4*(-t4+t3_aux);
continue_program=acept_solution(t0,t1_aux,t2_aux,t3_aux,t4,Vop,Vbk);

switch continue_program
    case 1
        [x_acceleration,x_cruising,x_coasting,x_breaking,d]=distance_calculator(t0,t1_aux,t2_aux,t3_aux,t4,  Vi, Vop, Vbk, Vf );
        [first_stop,second_stop,third_stop]=time_steps(t1_aux,t2_aux,t3_aux,time_vector);
        [v_vector_aux]=velocity_determination(Vi,t0,first_stop,second_stop,third_stop,t4,a1,a2,a3,a4,time_vector);
        
        for i=1:lengthTime_vector,
            if(i==1)
                aux=1;
                next_changekm=Vmax(aux+1,1)*km2m;
                actual_limit=Vmax(aux,2)*kmh2ms;
                next_limit=Vmax(aux+1,2)*kmh2ms;
                Vmax_vector(i,1)=actual_limit;
            elseif(i>1 && x_vector(i-1,1)<next_changekm)
                Vmax_vector(i,1)=actual_limit;
            elseif(x_vector(i-1,1)>next_changekm && (aux+1)<limits_number)
                aux=aux+1;
                next_changekm=Vmax(aux+1,1)*km2m;
                actual_limit=Vmax(aux,2)*kmh2ms;
                next_limit=Vmax(aux+1,2)*kmh2ms;
                Vmax_vector(i,1)=actual_limit;
            elseif(x_vector(i-1,1)>next_changekm && (aux+1)==limits_number)
                aux=aux+1;
                next_changekm=xf;
                actual_limit=Vmax(aux,2)*kmh2ms;
                next_limit=0;
                Vmax_vector(i,1)=actual_limit;
            end
            
            state_deb(i,1)=state;
            %TRAIN STATE MACHINE
            switch state
                case 0
                    v_vector(i,1)=Vi;
                    x_vector(i,1)=xi;
                    Vf=Vf;
                    state=1;
                    
                case 1  %ACCELERATION
                    %calculo da velocidade e posi��o em cada itera��o
                    v_vector(i,1)=(a1)*((time_vector(i,1)-time_vector(i-1,1)))+v_vector(i-1,1);
                    x_vector(i,1)=x_vector(i-1,1)+v_vector(i-1,1)*((time_vector(i,1)-time_vector(i-1,1)));
                    
                    if(Vop<=actual_limit)
                        disp(['Vop dentro do limite actual ', num2str(Vop), ' ', num2str(actual_limit)])
                        if(Vop<=next_limit && time_vector(i,1)>t1_aux)
                            state=2;
                        elseif(Vop>next_limit && time_vector(i,1)>t1_aux)
                            state=2;
                        elseif(Vop>next_limit)
                            disp('AQUI')
                            Vf=next_limit;
                            coasting_duration=(Vf^2-v_vector(i,1)^2)/(2*a3);
                            coasting_point=next_changekm-coasting_duration;
                            breaking_duration=(Vf^2-v_vector(i,1)^2)/(2*a4);
                            breaking_point=next_changekm-breaking_duration;
                            if(x_vector(i,1)>breaking_point && x_vector(i,1)<coasting_point && time_vector(i,1)>t2_aux)
                                state=4;
                            elseif(x_vector(i,1)>=coasting_point && time_vector(i,1)<t2_aux)
                                state=3;
                            end
                        end
                        
                    elseif(Vop>actual_limit)
                        disp(['Vop fora do limite actual ', num2str(Vop), ' ', num2str(actual_limit)])
                        if(abs(actual_limit-v_vector(i,1))>-V_tolerance && abs(actual_limit-v_vector(i,1))<V_tolerance)
                            state=2;
                        end
                    end
                    
                    
                case 2  %CRUISING
                    %calculo da velocidade e posi��o em cada itera��o
                    v_vector(i,1)=(a2)*((time_vector(i,1)-time_vector(i-1,1)))+v_vector(i-1,1);
                    x_vector(i,1)=x_vector(i-1,1)+v_vector(i-1,1)*((time_vector(i,1)-time_vector(i-1,1)));
                    
                    if(abs(actual_limit-v_vector(i,1))>-V_tolerance && abs(actual_limit-v_vector(i,1))<V_tolerance)
                        disp('CRUISING EM REGIME DE LIMITE DE VELOCIDADE')
                        if(actual_limit<next_limit)
                            %Aqui acelera
                            disp('est� a entrar aqui')
                            Vop=VopGen(actual_limit,next_limit);
                            Vi=v_vector(i,1);
                            Vf=0;
                            t0=time_vector(i,1);
                            [t1_aux,t2_aux,t3_aux]=times_calculation(Vi,Vf,Vop,t0,t4,a1,a2,a3,a4,Travel_d-x_vector(i,1));
                            Vbk=a4*(-t4+t3_aux);
                            while(acept_solution(t0,t1_aux,t2_aux,t3_aux,t4,Vop,Vbk)==0)
                                Vop=VopGen(actual_limit,next_limit);
                                Vi=v_vector(i,1);
                                Vf=0;
                                t0=time_vector(i,1);
                                tf=t4;
                                distancetoend=Travel_d-x_vector(i,1);
                                [t1_aux,t2_aux,t3_aux]=times_calculation(Vi,Vf,Vop,t0,tf,a1,a2,a3,a4,distancetoend);
                                Vbk=a4*(-t4+t3_aux);
                                [x_acceleration,x_cruising,x_coasting,x_breaking,d]=distance_calculator(t0,t1_aux,t2_aux,t3_aux,t4, Vi, Vop, Vbk, Vf);
                            end
                            if(x_vector(i,1)>next_changekm)
                                disp('GEROU NOVA VELOCIDADE E MANDOU ACELERAR')
                                state=1;
                            end
                                                        
                        elseif(actual_limit>next_limit)
                            disp('� necess�rio travar apos este cruising')
                            Vf=next_limit;
                            coasting_duration=(Vf^2-v_vector(i,1)^2)/(2*a3);
                            coasting_point=next_changekm-coasting_duration;
                            breaking_duration=(Vf^2-v_vector(i,1)^2)/(2*a4);
                            breaking_point=next_changekm-breaking_duration;
                            if(x_vector(i,1)>breaking_point && x_vector(i,1)<coasting_point && time_vector(i,1)>t3_aux)
                                state=4;
                            elseif(x_vector(i,1)>=coasting_point && time_vector(i,1)<t3_aux)
                                state=3;
                            end
                        end
                    elseif(v_vector(i,1)<actual_limit)
                        disp('trol')
                        if(time_vector(i,1)>t2_aux && actual_limit<=next_limit)
                            disp('depois do trol veio para o estado normal')
                            state=3;
                        elseif(actual_limit>next_limit && Vop>next_limit && v_vector(i,1)>next_limit)
                            disp('trol 2')
                            Vf=next_limit;
                            coasting_duration=(Vf^2-v_vector(i,1)^2)/(2*a3);
                            coasting_point=next_changekm-coasting_duration;
                            breaking_duration=(Vf^2-v_vector(i,1)^2)/(2*a4);
                            breaking_point=next_changekm-breaking_duration;
                            if(x_vector(i,1)>breaking_point && x_vector(i,1)>coasting_point && time_vector(i,1)>t2_aux)
                                state=4;  
                                disp('trol2 mandou travar')
                            elseif(x_vector(i,1)<coasting_point && x_vector(i,1)<breaking_point && time_vector(i,1)<t2_aux)
                                state=3;
                                disp('trol2 mandou coasting')
                            end
                        elseif(actual_limit>next_limit && Vop<next_limit)
%                             if(time_vector(i,1)>t2_aux)
%                                 state=3;
%                             end
                        end
                    end
                    
                    
                case 3  %COASTING
                    %calculo da velocidade e posi��o em cada itera��o
                    v_vector(i,1)=(a3)*((time_vector(i,1)-time_vector(i-1,1)))+v_vector(i-1,1);
                    x_vector(i,1)=x_vector(i-1,1)+v_vector(i-1,1)*((time_vector(i,1)-time_vector(i-1,1)));
  
                    disp('entrou em coasting')   

                    if(v_vector(i,1)<actual_limit && v_vector(i,1)<actual_limit)
                        if(time_vector(i,1)>t3_aux)
                            disp('entra aqui ##########################################')
                            %state=4;
                            Vf=next_limit;
                            coasting_duration=(Vf^2-v_vector(i,1)^2)/(2*a3);
                            coasting_point=next_changekm-coasting_duration;
                            breaking_duration=(Vf^2-v_vector(i,1)^2)/(2*a4);
                            breaking_point=next_changekm-breaking_duration;
                            if(x_vector(i,1)>breaking_point && x_vector(i,1)>coasting_point && time_vector(i,1)>t2_aux)
                                state=4;  
                                disp('trol2 mandou travar')
                                Vf
                            elseif(x_vector(i,1)<coasting_point && x_vector(i,1)<breaking_point && time_vector(i,1)<t2_aux)
                                state=3;
                                disp('trol2 mandou coasting')
                                Vf
                             end
                            
                        end
                        
                    elseif(v_vector(i,1)>actual_limit)
                        disp('entra aqui?????????????????????????????????????????????????')
                    end
                    
                    
                case 4  %BREAKING
                    %calculo da velocidade e posi��o em cada itera��o
                    v_vector(i,1)=(a4)*((time_vector(i,1)-time_vector(i-1,1)))+v_vector(i-1,1);
                    x_vector(i,1)=x_vector(i-1,1)+v_vector(i-1,1)*((time_vector(i,1)-time_vector(i-1,1)));
                    
                    if(Vf==0)
                        if(abs(v_vector(i,1)-0)>-V_tolerance && abs(v_vector(i,1)-0)<V_tolerance)
                            disp('ACABOU')
                        end
                    elseif(Vf>0)
                        if(abs(v_vector(i,1)-Vf)>-V_tolerance && abs(v_vector(i,1)-Vf)<V_tolerance)
                            state=2;
                        end
                    end
                    
                    if(v_vector(i,1)<0)
                        disp('estava a dar negativa')
                        v_vector(i,1)=0;
                    end
                    
                otherwise
                    disp('n�o contemplado')
                    
            end
            if(i==lengthTime_vector && abs(v_vector(lengthTime_vector,1)-0)>-V_tolerance && abs(v_vector(lengthTime_vector,1)-0)<V_tolerance && abs(xf-x_vector(lengthTime_vector,1))>-D_tolerance && abs(xf-x_vector(lengthTime_vector,1))<D_tolerance)
                disp('Percurso Encontrado')
            else
                while(acept_solution(t0,t1_aux,t2_aux,t3_aux,t4,Vop,Vbk)==0)
                    Vop=VopGen(0,v_max);
                    Vi=0*kmh2ms; %velocidade inicial
                    Vf=0*kmh2ms;
                    t0=0;
                    tf=t4;
                    [t1_aux,t2_aux,t3_aux]=times_calculation(Vi,Vf,Vop,t0,t4,a1,a2,a3,a4,Travel_d);
                    Vbk=a4*(-t4+t3_aux);
                    [x_acceleration,x_cruising,x_coasting,x_breaking,d]=distance_calculator(t0,t1_aux,t2_aux,t3_aux,t4, Vi, Vop, Vbk, Vf);
                end
                i=1;
                estado=0;
                
            end
            
            
        end
        
    otherwise
               
end

disp('Plotting...')

figure()
    subplot(2,1,1);
        plot(time_vector(:,1), v_vector(:,1)./kmh2ms, 'MarkerEdgeColor','b');
        grid on;
        legend('v_{train}');
        xlabel('m')
        ylabel('km/h')
        title('Velocity vs. Distance')
    subplot(2,1,2);
        plot(time_vector(:,1), x_vector(:,1));
        %fill( time_vector,  v_vector, 'blue');
        grid on;
        %axis([0 1 -15 0])
        legend('v_{train}');
        xlabel('s')
        ylabel('m/s')
        title('Velocity vs. Time')

        
figure()
    plot(x_vector(:,1), v_vector(:,1)./kmh2ms, x_vector(:,1), Vmax_vector(:,1)./kmh2ms);
    grid on;
    hold on;
    xlabel('m');
    ylabel('v/m');
    title('Distance')
    %legendCell(run)=cellstr(num2str(Vop_log(1,run)', 'V_{op}=%f\n'));
    
figure()
    subplot(2,1,1);
        plot(x_vector(:,1), v_vector(:,1)./kmh2ms, x_vector(:,1), v_vector_aux(:,1)./kmh2ms);
        grid on;
        legend('v_{train}');
        xlabel('m')
        ylabel('km/h')
        legend('state_machine', 'Previous')
        title('Velocity vs. Distance')
    subplot(2,1,2);
        plot(time_vector(:,1), v_vector(:,1)./kmh2ms, time_vector(:,1), v_vector_aux(:,1)./kmh2ms);
        %fill( time_vector,  v_vector, 'blue');
        grid on;
        %axis([0 1 -15 0])
        legend('v_{train}');
        xlabel('s')
        ylabel('m/s')
        legend('state_machine', 'Previous')
        title('Velocity vs. Time'); 
   

figure()
    subplot(3,1,1);
        plot(time_vector(:,1), state_deb(:,1));
        grid on;
        legend('v_{train}');
        xlabel('m')
        ylabel('km/h')
        legend('state_machine', 'Previous')
        title('Velocity vs. Distance')
    subplot(3,1,2);
        plot(time_vector(:,1), v_vector(:,1)./kmh2ms);
        %fill( time_vector,  v_vector, 'blue');
        grid on;
        %axis([0 1 -15 0])
        legend('v_{train}');
        xlabel('s')
        ylabel('m/s')
        legend('state_machine', 'Previous')
        title('Velocity vs. Time'); 
    subplot(3,1,3);
        plot(time_vector(:,1), x_vector(:,1));
        %fill( time_vector,  v_vector, 'blue');
        grid on;
        %axis([0 1 -15 0])
        legend('v_{train}');
        xlabel('s')
        ylabel('m/s')
        legend('state_machine', 'Previous')
        title('Velocity vs. Time');         
      
        
disp('END')