
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

% Loop through each file and display its name

load( files{1} )
G1 = oe(data , [8, 8, 1] ) ; 
Gf1 = spa(data, 8191, freqs) ; 

G = stack( 1 , G1 ) ; 
Gf = stack( 1 , Gf1 ) ; 

for idx = 2:length(files)
    disp(files{idx});
    load( files{idx} )
    G1 = oe(data , [8, 8, 1] ) ; 
    Gf1 = spa( data , 8191 , freqs ) ; 
    G = stack( 1 , G , G1 ) ; 
    Gf = stack( 1 , Gf , Gf1 ) ; 
end

%%



figure();
nyquistplot(Gf1, G, freqs, opts, 'sd ', 2.45);