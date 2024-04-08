clear all ; 
close all ; 
clc ; 

%% our systems transfert function 

G = tf( [ 4 , 2 ] , [ 1 , 4 , 4 ] ) ; 
mag = 5 ; 
R = tf( mag , [ 1 , 0 ] ) ; 






