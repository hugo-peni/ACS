
clear all ; 
close all ; 
clc ; 

load coeff.mat
%%


file_name = 'G.mat' ; 
current_folder_path = pwd ; 
full_file_path = fullfile( current_folder_path , file_name ) ; 

if( exist( full_file_path , 'file') == 2)
    load( 'G.mat' )
else
    load('Data/logs_7.mat')
    G = oe(data , [8, 8, 1] ) ;
    save('G.mat', 'G');
end

[B1, A1] = tfdata(G, 'v') ; 
T1s = G.Ts ; 

Ts=G.Ts; B=G.b; A=G.f;

disp( A1 )

%%

%Pideal = [0.8, 0.9 , 0.95] ; 
%Pideal = [1 , 0.95, 0.95 , 0.99] ; 
Pideal = [0.95, 0.95 , 0.99] ; 

Pchar = poly( Pideal ) ; 

Hs = [1 , -1 ] ; 
Hr = [1 , 1 ] ; 

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

T = sum( Pchar ) / sum( B )

CL=tf(conv(T,B),Pplaced,Ts,'variable','z^-1')

[Y,time] = step(CL) ;
crop = length( time ) / 20 ; 

figure(1)
%plot( time(1:crop) , Y(1:crop) )
step( CL )

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













