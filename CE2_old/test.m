clear;
close;
clc;

load("data1.mat");


G1 = bj(data1, [5, 5, 5, 5, 1]);
Gf1 = spa(data1, 400);

figure(1)
bode(Gf1)
title('Bode diagram of the first system')

%% Plot Nyquist 
P = nyquistoptions; 
P.ConfidenceRegionDisplaySpacing=1; 
P.ShowFullContour='off';
figure();
subplot(121) 
nyquistplot(G1,Gf1,P,'sd',2); 
title('Nyquist Diagram of the first system') ; 
legend('From parmatric model', 'From non-parametric model');
subplot(122)
P.Xlim=[45, 50]; P.Ylim=[-30, 5]; 
h=nyquistplot(G1,Gf1,P,'sd',2); 
title('Zoom in Nyquist Diagram') ; 
axis equal
