clear ; close ; clc ;

%compute W2
load('data1.mat') ; load('data2.mat') ; load('data3.mat') ; load('data4.mat') ; 
G1=bj(data1, [5,5,5,5,1]) ; G2=bj(data2, [5,5,5,5,1]) ; G3=bj(data3, [5,5,5,5,1]) ;G4=bj(data4, [5,5,5,5,1]);
Ts = data1.Ts ;

%load K_red and convert to tf
load('result_2_3.mat') ; 
[num, den] = ss2tf(K_red.A,K_red.B,K_red.C,K_red.D) ; 
% K_red = c2d(tf(num,den),Ts) ; 

G = cat(1,G1,G2,G3,G4) ; Gmm = stack(1, G1, G2, G3, G4);

% Frequency grid
W = logspace(-2,log10(pi/Ts),1000); 

% Define the weighting filter W1(s) in continuous time
s = tf('s') ; m = 0.5 ; omega_b = 370 ; 
W1s = m * (s+omega_b) / (s+1e-5) ; 
% Convert to discrete time using zero-order hold
W1z = c2d(W1s, Ts, 'zoh'); 

% W1z = 1/makeweight(0.01,3,2,Ts);

CON = struct('W1',W1z,'W2',W2,'W3', []); % H infinity constraints
OBJ = struct('W1',1,'W2',1 ,'W3',[],'norm','inf'); % Objectives to minimize
opt = sdpsettings('solver','mosek','verbose',0); % YALMIP settings

PAR = datadrivenOptions('Kc',K_red,'order',5,...
                        'W',W, ...
                        'force_integrator',0,... % if 'force_integrator',1, best practice is to already have an integrator in Kc
                        'sdpsettings',opt,...
                        'max_iteration',10); 

%%                    
K_optimal = datadriven_ACS(Gmm, OBJ, CON, PAR) ; 

%%

T_optimal = feedback(G1*K_optimal,1) ; 
step(T_optimal) ; 

