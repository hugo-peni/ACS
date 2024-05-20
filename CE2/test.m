
clear all ; 
close all ; 
clc ; 

%% 

Ts = 0.002 ; % sampling time in s 

s=tf('s'); m=0.505; omega_b=2;
W1s=m*(s+omega_b)/(s+1e-5) ;
W1d=c2d(W1s,Ts,'tustin') ;

figure(1)
bode( 1 / W1s )

%%