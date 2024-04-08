clear all ; 
close all ; 
clc ; 

%% define the transfert functions for system 1 

G = tf( [2 , -1 ] , [ 1 , 0 , 2 ]) ;
K = tf( 1 , 1 ) ; 

G1 = minreal( K * G / ( 1 + K * G )) ; 

P = pole( G1 )
stability = isstable( G1 )

%% define transfert functions for system 2 

K2 = tf( [ 4 , 2 ] , [ 2 , -1 ] ) ; 

G2 = minreal( ( K2 * G) / ( 1 + ( K2 * G ) ) )
P = pole( G2 )
stability2 = isstable( G2 )

%% 

G3 = tf( [ 4 , 2 ] , [ 1 , 4 , 4 ] )

figure(1)
bode( G3 )
hold on 
bode( G2)
hold off
legend( "C_G2" , "G3" )

%%
figure(1)
step( G1 )

figure(2)
step( G3 )