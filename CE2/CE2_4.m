
clc
K_pid_init = pid(1e-2, 1e-5 ,1e-5, Ts, Ts);

%fprintf('Settling times with initial pid controller : %.2fs, %.2fs, %.2fs, %.2fs \n',stepinfo(feedback(G1*K_pid_init,1)).SettlingTime,stepinfo(feedback(G2*K_pid_init,1)).SettlingTime,stepinfo(feedback(G3*K_pid_init,1)).SettlingTime,stepinfo(feedback(G4*K_pid_init,1)).SettlingTime)
% step(feedback(Gmm*K_pid_init,1)) ; grid on ; 

%%

CON_pid = struct('W1',W1s,'W2',W3,'W3', []) ; % H infinity constraints
OBJ_pid = struct('W1',1,'W2',[] ,'W3',[],'norm','two') ; % Objectives to minimize


opt = sdpsettings('solver','mosek','verbose',0) ; % YALMIP settings
PAR_pid = datadrivenOptions( 'Kc' , K_pid_init , 'order' , 2 , ...
                        'W' , W , ... 
                        'force_integrator' , 0 ,... 
                        'sdpsettings' , opt , ...
                        'max_iteration' , 10 ) ; 


%K_pid = datadriven_ACS( G ,OBJ_pid , CON_pid , PAR_pid ) ;
