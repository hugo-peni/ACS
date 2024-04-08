clear ; close ; clc ;

% Load data & models
load('data1.mat') ; load('data2.mat') ; load('data3.mat') ; load('data4.mat') ; 
G1=bj(data1, [5,5,5,5,1]) ; G2=bj(data2, [5,5,5,5,1]) ; G3=bj(data3, [5,5,5,5,1]) ;G4=bj(data4, [5,5,5,5,1]);
Gmm = stack(3, G1, G2, G3, G4) ; 

% Compute W2 in discrete time
Gnom = G1 ; 
[Gu, info] = ucover(Gmm, Gnom, 4) ;
W2=info.W1 ; 
% Define the weighting filter W1(s) in continuous time
s = tf('s') ; m = 0.5 ; omega_b = 10 ; 
W1s = m * (s+omega_b) / (s+1e-5) ; 
% Convert to discrete time using zero-order hold
Ts = 0.01 ; 
W1z = c2d(W1s, Ts, 'zoh'); 

% Check the magnitude of W1^-1(s) at high frequencies
[mag, ~] = bode(W1s^-1, 100); magdB = 20*log10(mag) ; 
% fprintf('Magnitude of W1(s)^-1 at high frequencies : %.2f dB\n',magdB) ; 

% Reduce magnitude of u(t) with W3
W3 = 4.5 ; 

%Compute controller K 
K = mixsyn(tf(Gnom),W1z,W3,W2) ; 

% Define the CL transfer functions
S=feedback(1,Gmm*K) ; T=feedback(Gmm*K,1) ; U=feedback(K,Gmm) ; 

% Plot CL step response
% figure ; bodemag(S) ; grid on ; title('Sensitivity functions S')
% figure ; bodemag(U) ; grid on ; title('Input sensitivity functions U')
% figure ; step(T) ; grid on ; legend('Step response closed loop') ;
% figure ; step(U) ; grid on ; legend('Control signal') ; 

% figure;
% subplot(2, 2, 1);bodemag(S);grid on;title('Sensitivity functions S');
% subplot(2, 2, 2);bodemag(U);grid on;title('Input sensitivity functions U');
% subplot(2, 2, 3);step(T);grid on; legend('Step response closed loop');
% subplot(2, 2, 4);step(U);grid on;legend('Control signal');

fprintf('Settling times with initial controller : %.2fs, %.2fs, %.2fs, %.2fs \n',stepinfo(feedback(G1*K,1)).SettlingTime,stepinfo(feedback(G2*K,1)).SettlingTime,stepinfo(feedback(G3*K,1)).SettlingTime,stepinfo(feedback(G4*K,1)).SettlingTime)
fprintf('Overshoot with initial controller : %.2f, %.2f, %.2f, %.2f \n',stepinfo(feedback(G1*K,1)).Overshoot,stepinfo(feedback(G2*K,1)).Overshoot,stepinfo(feedback(G3*K,1)).Overshoot,stepinfo(feedback(G4*K,1)).Overshoot)
%%
%Reduce the order of the controller 
% pzmap(K)
K_red = reduce(K,5) ; % reduce to order 5

% CL transfer functions with reduced controller
S_red=feedback(1,Gmm*K_red) ; T_red=feedback(Gmm*K_red,1) ; U_red=feedback(K_red,Gmm) ; 

% % figure ; step(T_red) ; grid on ; legend('Step response closed loop with reduced controller') ;
% figure ; step(U_red) ; grid on ; legend('Control signal with reduced controller') ; 
% figure ; step(T) ; grid on ; hold on ; step(T_red) ; legend('Initial controller','Reduced controller') ; 
% figure;
% subplot(1, 2, 1);step(U_red) ; grid on ; legend('Control signal with reduced controller') ; 
% subplot(1, 2, 2);step(T_red);grid on;legend('Step response with reduced controlelr');

fprintf('Settling times with controller %d : %.2fs, %.2fs, %.2fs, %.2fs \n',order(K_red),stepinfo(feedback(G1*K_red,1)).SettlingTime,stepinfo(feedback(G2*K_red,1)).SettlingTime,stepinfo(feedback(G3*K_red,1)).SettlingTime,stepinfo(feedback(G4*K_red,1)).SettlingTime)
fprintf('Overshoot with controller %d : %.2f, %.2f, %.2f, %.2f \n',order(K_red),stepinfo(feedback(G1*K_red,1)).Overshoot,stepinfo(feedback(G2*K_red,1)).Overshoot,stepinfo(feedback(G3*K_red,1)).Overshoot,stepinfo(feedback(G4*K_red,1)).Overshoot)
%% Closed loop norm of the objective 
%Compute the CL norm
cl_sys_nom = [W1z * feedback(1,Gnom*K_red) ; W3*feedback(K_red,Gnom)] ; 
cl_sys = [W1z * S_red ; W3*U_red] ; 

cl_norm_nom = norm(cl_sys_nom, Inf) ; 
fprintf('Closed loop norm with nominal model : %.4f\n', cl_norm_nom) ; 
cl_norm_mm = norm(cl_sys, Inf) ;
fprintf('Closed loop norm with multimodel : %.4f\n', cl_norm_mm) ;
%% Ex 2.4
G = cat(3,G1,G2,G3,G4) ; 

% Frequency grid
W = logspace(-2,log10(pi/Ts),1000); 

CON = struct('W1',[],'W2',[],'W3', []); % H infinity constraints
OBJ = struct('W1',W1s,'W2',W3 ,'W3',[],'norm','inf'); % Objectives to minimize
opt = sdpsettings('solver','mosek','verbose',0); % YALMIP settings

PAR = datadrivenOptions('Kc',K_red,'order',5,...
                        'W',W, ...
                        'force_integrator',0,... 
                        'sdpsettings',opt,...
                        'max_iteration',10); 
                    
K_dd = datadriven_ACS(G, OBJ, CON, PAR) ; 
%%
S_dd=feedback(1,Gmm*K_dd) ; T_dd=feedback(Gmm*K_dd,1) ; U_dd=feedback(K_dd,Gmm) ; 
figure ; step(U_dd) ; grid on ; legend('Control signal with datadriven controller') ; 
figure ; step(T) ; grid on ; hold on ; step(T_red) ; step(T_dd) ; legend('Initial controller','Reduced controller','Datadriven controller') ; 
fprintf('Settling times with reduced controller : %.2fs, %.2fs, %.2fs, %.2fs \n',stepinfo(feedback(G1*K_dd,1)).SettlingTime,stepinfo(feedback(G2*K_dd,1)).SettlingTime,stepinfo(feedback(G3*K_dd,1)).SettlingTime,stepinfo(feedback(G4*K_dd,1)).SettlingTime)
fprintf('Overshoot with datadriven controller : %.2f, %.2f, %.2f, %.2f \n',stepinfo(feedback(G1*K_dd,1)).Overshoot,stepinfo(feedback(G2*K_dd,1)).Overshoot,stepinfo(feedback(G3*K_dd,1)).Overshoot,stepinfo(feedback(G4*K_dd,1)).Overshoot)
%% 2.4.2
clc
K_pid_init = pid(1e-2, 1e-5 ,1e-5, Ts, Ts);

fprintf('Settling times with initial pid controller : %.2fs, %.2fs, %.2fs, %.2fs \n',stepinfo(feedback(G1*K_pid_init,1)).SettlingTime,stepinfo(feedback(G2*K_pid_init,1)).SettlingTime,stepinfo(feedback(G3*K_pid_init,1)).SettlingTime,stepinfo(feedback(G4*K_pid_init,1)).SettlingTime)
% step(feedback(Gmm*K_pid_init,1)) ; grid on ; 
%%
CON_pid = struct('W1',W1s,'W2',W3,'W3', []); % H infinity constraints
OBJ_pid = struct('W1',1,'W2',[] ,'W3',[],'norm','two'); % Objectives to minimize

opt = sdpsettings('solver','mosek','verbose',0); % YALMIP settings
PAR_pid = datadrivenOptions('Kc',K_pid_init,'order',2,...
                        'W',W, ...
                        'force_integrator',0,... 
                        'sdpsettings',opt,...
                        'max_iteration',10); 
                    
K_pid = datadriven_ACS(Gmm,OBJ_pid,CON_pid,PAR_pid);
%%
S_pid=feedback(1,Gmm*K_pid) ; T_pid=feedback(Gmm*K_pid,1) ; U_pid=feedback(K_pid,Gmm) ; 
figure ; step(U_pid) ; grid on ; legend('Control signal with datadriven controller') ; 
figure ; step(T) ; grid on ; hold on ; step(T_pid) ; legend('Initial controller','DD PID controller') ; 
fprintf('Settling times with dd pid controller : %.2fs, %.2fs, %.2fs, %.2fs \n',stepinfo(feedback(G1*K_pid,1)).SettlingTime,stepinfo(feedback(G2*K_pid,1)).SettlingTime,stepinfo(feedback(G3*K_pid,1)).SettlingTime,stepinfo(feedback(G4*K_pid,1)).SettlingTime)
fprintf('Overshoot with dd pid controller : %.2f, %.2f, %.2f, %.2f \n',stepinfo(feedback(G1*K_pid,1)).Overshoot,stepinfo(feedback(G2*K_pid,1)).Overshoot,stepinfo(feedback(G3*K_pid,1)).Overshoot,stepinfo(feedback(G4*K_pid,1)).Overshoot)

