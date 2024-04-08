clear ; close ; clc ; 

%load data
load("data1.mat")
G1 = bj(data1, [5, 5, 5, 5, 1]);
G1_C = d2c(G1);
[A,B,C,D] = ssdata(G1_C);
%% LMI
% Define the decision variables
n = size(A,1);
L = sdpvar(n, n,'symmetric'); % symmetric positive definite matrix
Gamma = sdpvar(1,1) ; 
Y = sdpvar(1, n);

% Define the objective function
objective = trace(C*L*C') + trace(Gamma) ;  

% Define the LMI constraints
LMI1 = [A*L + L*A' - B*Y - Y'*B' + B*B']<= 0;
LMI2 = [Gamma Y ; Y' L] >= 0;

% Solve the optimization problem
constraints = [LMI1, LMI2,L >= 0];
options = sdpsettings('solver', 'sedumi','Verbose',0);
sol = optimize(constraints, objective, options);
Y = value(Y) ; L = value(L) ; K_lmi = Y /L 
%% LQR 
% Define the weighting matrices
Q = C'*C;
R = 1;
% Compute the LQR controller gain
K_lqr = lqr(A, B, Q, R)
err_k = norm(K_lmi-K_lqr,2) 
%% Plots
%step response of closed loop system
sys_lmi = ss(A - B*K_lmi,B,C,D) ;
step(sys_lmi) ; grid on ; 
% exportgraphics(gcf,'Figures/stepCL.png','Resolution',300) ; 

