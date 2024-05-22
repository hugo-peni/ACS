
clear all ; 
close all ; 
clc ; 


files = {'Data/logs_1.mat', 'Data/logs_3.mat', 'Data/logs_5.mat', ...
         'Data/logs_7.mat', 'Data/logs_9.mat', 'Data/logs_11.mat'};


current_folder_path = pwd ; 
full_file_path = fullfile( current_folder_path , 'G.mat' ) ; 

if( exist( full_file_path , 'file') == 0 )
    load( files{1} )
    G0 = oe(data , [8, 8, 1] ) ;
    G = stack( 1 , G0 ) ;
    for idx = 2:length(files)
        load( files{idx} )
        G0 = oe(data , [8, 8, 1] ) ;
        G = stack( 1 , G , G0 ) ;
    end
    save( 'G.mat' , 'G' );
end

clear current_folder_path
clear files 
load( 'G.mat' )

chosen_nominal_model = 6 ; 
nom_model = G(:,:,chosen_nominal_model,1) ; % take model 4 as the nominal model 

[B, A] = tfdata( nom_model , 'v' ) ; 
Ts = nom_model.Ts ; 


%%

type = "poles_1" ; 

if( type == "poles_1" )
    Pideal = poly([0.8, 0.9 , 0.95]) ;
elseif( type == "poles_2")
    Pideal = poly([0.95, 0.95 , 0.99]) ; 
end

Hs = [1 , -1 ] ; 
%Hr = 1 ; 
Hr = [1 , 1 ] ; 


[R,S]=poleplace(B,A,Hr,Hs, Pideal ) ; 
Pplaced = conv(A,S)+conv(B,R) ; 

disp( Pideal )
disp( Pplaced(1:4) )


T = sum( Pplaced ) / sum( B ) ; 

CL = tf(conv(T,B),Pplaced,Ts,'variable','z^-1') ; 

figure(1)
step( CL )

U = tf(conv(A,R),Pplaced,Ts,'variable','z^-1') ; 

figure(2)
step( U )


%% Q-parametrization to complete 

s = tf('s') ; m = 0.5 ; omega_b = 10 ; 
W1s = m * (s+omega_b) / (s+1e-5) ; 
% Convert to discrete time using zero-order hold
Ts = 0.01 ; 
W1z = c2d(W1s, Ts, 'zoh'); 

% Compute W2 in discrete time
[~, info] = ucover(G, G(:,:,chosen_nominal_model,:) , 7) ;
W2=info.W1 ; 


fprintf('\n') ; 
orderQinit = 9 ; %Can be changed (between 3 and 6 maybe) 
Q_init = -1 + 2.*rand(1,orderQinit+1) ; 


%%
options = optimoptions('fmincon', 'Display', 'iter'); 
options.MaxFunctionEvaluations = Inf ; 
options.MaxIterations = 10000 ; 
obj = @(Q) define_objective( A , B , S , Q , Pplaced , W1z , W2 , Ts ) ;

%call solver
[Q,~,~,OUTPUT] = fmincon(obj, Q_init, [], [], [], [], [], [], @constraint, options);
Q_best = OUTPUT.bestfeasible

 

%% 


function obj = define_objective(A,B,S,Q,P,W1,W2,Ts)
    num = conv(A,diff_coeff(S,conv(B,Q))) ; 
    den = P ; 
    SensS = tf(num,den,Ts, 'variable','z^-1') ;
    obj = norm(W1*SensS,Inf) + norm(W2*(1-SensS),Inf) ; 
end

function SensS_Inf = define_SensS(A,B, S,Q,P, Ts)
    %returns the inf norm of SensS in dB
    num = conv(A,diff_coeff(S,conv(B,Q))) ; 
    den = P ; 
    SensS_Inf = 20*log10(norm(tf(num,den,Ts, 'variable','z^-1'),Inf)) ;
end

function SensU_Inf = define_SensU(A,R, Q,P, Ts)
    %return the inf norm of SensU in dB
    num = conv(A,diff_coeff(R,-conv(A,Q))) ; 
    den = P ; 
    SensU_Inf = 20*log10(norm(tf(num,den,Ts, 'variable','z^-1'),Inf)) ;
end

function diff = diff_coeff(array1, array2)
    % Determine the lengths of the input arrays
    len1 = length(array1); len2 = length(array2);
    % Calculate the difference in lengths
    diff_len =abs((len1 - len2));
    
    if diff_len == 0
        diff = array1-array2 ;
    % Pad the shorter array with zeros
    elseif len1 > len2
        diff  = array1 - [array2, zeros(1, diff_len)]; 
    else
        diff = [array1, zeros(1, diff_len)] - array2; 
    end
end


function [c, ceq] = constraint(Q, A, B, S, R, Ts ) 
    SensS_Inf = define_SensS(A, B, S , Q, P , Ts);
    SensU_Inf = define_SensU(A, R , Q, P , Ts);
    c = [SensS_Inf - 6; SensU_Inf + 10];  % Inequality constraints
    %S<6dB & U<-10dB
    ceq = [];                                               % No equality constraints
end


