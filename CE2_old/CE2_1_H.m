clear ; close ; clc ; 

%% Load Data

load("data1.mat")
load("data2.mat")
load("data3.mat")
load("data4.mat")


%% System identification 


G1 = bj(data1, [5, 5, 5, 5, 1]);
Gf1 = spa(data1, 400);
G2 = bj(data2, [5, 5, 5, 5, 1]);
Gf2 = spa(data2, 400);
G3 = bj(data3, [5, 5, 5, 5, 1]);
Gf3 = spa(data3, 400);
G4 = bj(data4, [5, 5, 5, 5, 1]);
Gf4 = spa(data4, 400);

%% Frequency response

figure(1)
bode(Gf1);

figure(2)
bode(Gf2);

figure(3)
bode(Gf3);

figure(4)
bode(Gf4);



%% Nyquist plots 

P = nyquistoptions;
P.ConfidenceRegionDisplaySpacing=1; 
P.ShowFullContour='off';
figure(5);
subplot(121) 
nyquistplot(Gf1,G1,P,'sd',2); 
subplot(122)
P.Xlim=[45, 50]; 
P.Ylim=[-30, 5]; 
h=nyquistplot(Gf1,G1,P,'sd',2); 
title('Zoom in Nyquist Diagram') 
axis equal


%% Wheighting filter that converts multimodal uncertainty into multiplicative one

Gmm = stack(1, G1, G2, G3, G4);

Gnom1 = G1;
[Gu1, info1] = ucover(Gmm, Gnom1, 1);
W2_1 = info1.W1opt;
E2_1 = Gmm/Gnom1 -1 ; 

Gnom2 = G2;
[Gu2, info2] = ucover(Gmm, Gnom2, 1);
W2_2 = info2.W1opt;
E2_2 = Gmm/Gnom2 -1 ; 
 
Gnom3 = G3;
[Gu3, info3] = ucover(Gmm, Gnom3, 1);
W2_3 = info3.W1opt;
E2_3 = Gmm/Gnom3 -1 ; 

Gnom4 = G4;
[Gu4, info4] = ucover(Gmm, Gnom4, 1);
W2_4 = info4.W1opt;
E2_4 = Gmm/Gnom4 -1 ; 

%% Bode plot of the filters 

figure(6)
bode(W2_1);
hold on
bode(W2_2);
bode(W2_3);
bode(W2_4);
hold off


figure(7)
bode(E2_1);
hold on
bode(E2_2);
bode(E2_3);
bode(E2_4);
hold off
