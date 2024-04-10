
clear all ; 
close all ; 
clc ; 

%% miscellaneous varibales and plotting setup 

Ts = 0.002 ; % sampling time in s 
freqs = (pi/4096:pi/4096:pi) / Ts ; 

opts = nyquistoptions; 
opts.ConfidenceRegionDisplaySpacing = 3;
opts.ShowFullContour = 'off';

%% load data and identify the model 

% Define the list of files
files = {'Data/logs_1.mat', 'Data/logs_3.mat', 'Data/logs_5.mat', ...
         'Data/logs_7.mat', 'Data/logs_9.mat', 'Data/logs_11.mat'};

% intitialization (identify model with first file)
load( files{1} )
G_current = oe(data , [8, 8, 1] ) ; 
Gf_current = spa(data, 8191, freqs) ; 
G = stack( 1 , G_current ) ; 
Gf = stack( 1 , Gf_current ) ; 

for idx = 2:length(files)
    disp(files{idx});
    load( files{idx} ) ; 
    G_current = oe(data , [8, 8, 1] ) ; 
    Gf_current = spa( data , 8191 , freqs ) ; 
    G = stack( 1 , G , G_current ) ; 
    Gf = stack( 1 , Gf , Gf_current ) ; 
end

%%

figure(1)

for i = 1:length(G)
    subplot( 2 , 3 , i )
    nyquistplot( G(:,:,i,:) , Gf(:,:,i,:) , opts , 'sd' , 2.45 ) ; 
end

%%

[~, info1] = ucover( G , G(:,:,1,:) , 4 ) ;
W_current = info1.W1 ; 
n_inf_current = norm( info1.W1 , inf ) ; 

W = stack( 1 , W_current ) ; 
n_inf = n_inf_current ; 

for i = 2:length(G)
    [~, info1] = ucover( G , G(:,:,i,:) , 4 ) ;
    W_current = info1.W1 ;  
    n_inf_current = norm( info1.W1opt , inf ) ; 
    W = stack( 1 , W , W_current ) ; 
    n_inf = [ n_inf , n_inf_current ] ;
end

%%

figure(2)
for i = 1:length(G)
    nom_S = strcat( 'G / Gnom' , num2str( i ) , ' - 1' ) ;
    W_S = strcat( 'W_2_' , num2str( i ) ) ;
    subplot(3, 2, i)
    bodemag( G / G(:,:,i,:) - 1 , W(:,:,i,:) ) ; grid on ; legend( nom_S , W_S ) ; 
    title( strcat( 'norm_inf(', W_S , ') = ', num2str( n_inf(i) ) ) )
end





