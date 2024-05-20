
clear all ; 
close all ; 
clc ; 

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

clc 

%Pideal = [0.8, 0.9 , 0.95] ; 
Pideal = [0.95, 0.95 , 0.99] ; 

Hs = [1, -1] ; 
Hr = 1 ;

[R,S]=poleplace(B,A,Hr,Hs,Pideal) ; 

fprintf("--R----\n")
disp(R)
fprintf("--S----\n")
disp(S)


% clc ; 
% 
Pplaced = conv(A,S)+conv(B,R) ; 
disp( Pplaced(1:3) )
% 
% disp( Pideal )
% disp( Pplaced(1:3) )
% 
% print("R = ")
% disp( R )
% print(" S = ")
% disp( S )





T = R(1)
T = Pplaced(1:3) / R(1) 

CL=tf(conv(T,B),Pplaced,Ts,'variable','z^-1')

[Y,time] = step(CL) ;
crop = length( time ) / 20 ; 

figure(1)
plot( time(1:crop) , Y(1:crop) )

figure(2)
step(CL)











