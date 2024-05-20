
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


% nom_model = G1 ; 
% 
% [B, A] = tfdata( nom_model , 'v' ) ; 
% T1s = nom_model.Ts ; 


%%

Pideal = [0.8, 0.9 , 0.95] ; 
%Pideal = [1 , 0.95, 0.95 , 0.99] ; 
%Pideal = [0.95, 0.95 , 0.99] ; 

Pchar = poly( Pideal ) ; 

Hs = [1 , -1 ] ; 
Hr = [1 , 1 ] ; 
%Hr = 1 ; 

% Hs = [ 1 , Beta , 1 ]
% Hr = [ 1 , Beta ] ; 

% [R,S]=poleplace(B,A,Hr,Hs,Pideal) ; 
[R,S]=poleplace(B,A,Hr,Hs, Pchar ) ; 

fprintf("--R----\n")
disp(R)
fprintf("--S----\n")
disp(S)


% clc ; 
% 
Pplaced = conv(A,S)+conv(B,R) ; 
disp( Pplaced(1:3) )

%T = R(1)
%T = Pideal / B(2)

%T = sum( Pchar ) / sum( B )
T = Pchar(1) / B(2)

CL=tf(conv(T,B),Pchar,Ts,'variable','z^-1')

[Y,time] = step(CL) ;
crop = length( time ) / 20 ; 

figure(1)
%plot( time(1:crop) , Y(1:crop) )
step( CL )

U=tf(conv(A,R),Pchar,Ts,'variable','z^-1')

figure(2)
step( U )


%%

Fs = 1 / Ts ; 
L = length( Y(1:crop) ) ; 

f = Fs/L*(0:(L/2)) ;

Y = fft(Y(1:crop)) ;

P2 = abs(Y/L) ;
P1 = P2(1:L/2+1) ;
P1(2:end-1) = 2*P1(2:end-1) ;

figure(90)
plot(f,P1,"LineWidth",3) 


f_max = f( find( P1 == max( P1 ) ) )

Beta = -2 * cos( 2*pi * f_max / Fs ) ; 

save("coeff.mat", "Beta")


%% Q-parametrization
% exportgraphics(gcf,'Figures/CompInitImpr.png','Resolution',500)

%% Q-parametrization

%compute the filters W1 and W2
% Define the weighting filter W1(s) in continuous time
s = tf('s') ; m = 0.5 ; omega_b = 10 ; 
W1s = m * (s+omega_b) / (s+1e-5) ; 
% Convert to discrete time using zero-order hold
Ts = 0.01 ; 
W1z = c2d(W1s, Ts, 'zoh'); 

% Load data & models
load('data1.mat') ; load('data3.mat') ; load('data4.mat') ; 
G1=bj(data1, [5,5,5,5,1]) ; G3=bj(data3, [5,5,5,5,1]) ;G4=bj(data4, [5,5,5,5,1]);
Gmm = stack(3, G1, G2, G3, G4) ; 
% Compute W2 in discrete time
[~, info] = ucover(Gmm, G2, 4) ;
W2=info.W1 ; 

%% Q-parametrization

fprintf('\n') ; 
orderQinit = 4 ; %Can be changed (between 3 and 6 maybe) 
Q_init = -1 + 2.*rand(1,orderQinit+1) ;  
options = optimoptions('fmincon', 'Display', 'iter'); 
options.MaxFunctionEvaluations = Inf ; 
options.MaxIterations = 10000 ; 
obj = @(Q) define_objective( A , B , S_improved , Q , P_improved , W1z , W2 , Ts ) ;
 


