%{
Data-driven command:
Data-driven mixed-sensitivity controller design. SISO only.
Computes the optimal controller for the mixed synthesis problem.
For the feedback loop

                    e       u
           r --->O--->[ K ]-->[ G ]---+---> y
               - |                    |
                 +<-------------------+

   the mixed-sensitivity design seeks a controller K that minimizes

                  || OBJ.W1*S ||
            min   || OBJ.W2*U ||
                  || OBJ.W3*T ||_norm
       subject to
       ||  CON.W1*S  ||_\infty ≤ 1
       ||  CON.W2*U  ||_\infty ≤ 1
       ||  CON.W3*T  ||_\infty ≤ 1

   where K = X/Y and
     * S = 1/(I+G*K) is the sensitivity function
     * U is the transfer from r to u (control effort)
     * T = I-S = G*K/(I+G*K) is the complementary sensitivity
     * norm is the norm specified in options (either 'two' or 'inf')
If multi-model uncertainity is used, use G = cat(1,G1,G2,..);


In the example below, compute optimal controller (of order 2)  solution to 

            min   || OBJ.W2*U ||_inf
       subject to
       ||  CON.W1*S  ||_\infty ≤ 1


%}

load data1
load data2
load data3
load data4;

Ts = data1.Ts;
Gf1 = spa(data1, 400);
G1 = bj(data1, [5, 5, 5, 5, 1]);
G2 = bj(data2, [5, 5, 5, 5, 1]);
G3 = bj(data3, [5, 5, 5, 5, 1]);
G4 = bj(data4, [5, 5, 5, 5, 1]);

Gmm = stack(1, G1, G2, G3, G4);

%% Set-up problem

% load mode

W1 = 1/makeweight(0.01,3,2,Ts);
CON = struct('W1',W1,'W2',W3,'W3', []); % H infinity constraints

OBJ = struct('W1',1,'W2',[] ,'W3',[],'norm','two'); % Objectives to minimize


% Some additional parameters
%Kc = tf(0.001,1,Ts); % initial stabilizing controller
Kc = pid(0, -1, 0, 0, Ts);
W = logspace(-2,log10(pi/Ts),1000); % Frequency grid at which we solve the problem

opt = sdpsettings('solver','mosek','verbose',0); % YALMIP settings

PAR = datadrivenOptions('Kc',Kc,'order',3,...
                        'W',W, ...
                        'force_integrator',0,... % if 'force_integrator',1, best practice is to already have an integrator in Kc
                        'sdpsettings',opt,...
                        'max_iteration',10); 
                    
% datadrivenOptions(name1,value1, name2, value2, ...)
% Structure of parameters used in datadriven command.
% --------------------------------------------------------------------
%     name          ->          description
% --------------------------------------------------------------------
% 'K'               -> intitial stabilizing controller
% 'order'           -> order of the final controller
% 'force_integrator'-> forces an integrator in the final controller
% 'W'               -> Frequency grid where the problem is solved
% 'max_iteration'   -> number of iteration before stopping
% 'sdpsettings'     -> sdpsettings structure for YALMIP
% --------------------------------------------------------------------
%
% Example use:
% options = datadrivenOptions('Kc',Kinit, 'order', orderK, ...
%                             'W',W, ...
%                             'force_integrator',true, ...
%                             'max_iteration',20);



%% Compute optimal controller

K_datadriven = datadriven_ACS(Gmm,OBJ,CON,PAR);


S_dd = feedback(1,Gmm*K_datadriven);

figure(1)
bodemag(S_dd,1/W1,W);
legend('Sensitivity','Constraints on sensitivity')

figure(2)
step(S_dd) 
grid on 
legend('Control signal with datadriven PID controller') 

shg
