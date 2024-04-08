clear ; close ; clc ; 

load("data1.mat")

G1 = bj(data1, [5, 5, 5, 5, 1]);
G1_C = d2c(G1);

[A,B,C,D] = ssdata(G1_C);

% Define the decision variables
n = size(A,1);
K = sdpvar(n, n); % state feedback gain matrix
L = sdpvar(n, n); % symmetric positive definite matrix

% Define the objective function
y1 = C*sdpvar(n,1);
y2 = -K*sdpvar(n,1);
Y = [y1 ; y2] ; 
objective = norm(Y,2) ; 

% Define the LMI constraints
LMI1 = [L C*A+B'*K; (C*A+B'*K)' -L] <= 0;
LMI2 = L >= eye(n);

% Solve the optimization problem
constraints = [LMI1, LMI2];
options = sdpsettings('solver', 'mosek');
sol = optimize(constraints, objective, options);

% Check the solution status
if sol.problem ~= 0
    error('Error in optimization: %s', sol.info)
end

% Retrieve the optimal solution
K_opt = value(K);
S_opt = value(S);