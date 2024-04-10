
clear ; close ; clc ;

%%

load('data1.mat')
G1 = bj(data1, [5, 5, 5, 5, 1]) ;  %parametric
Gf1 = spa(data1, 400) ;             %non parametric

%% Plot Nyquist 
% P = nyquistoptions; 
% P.ConfidenceRegionDisplaySpacing=1; 
% P.ShowFullContour='off';
% figure();
% subplot(121) 
% nyquistplot(G1,Gf1,P,'sd',2); 
% subplot(122)
% P.Xlim=[45, 50]; P.Ylim=[-30, 5]; 
% h=nyquistplot(G1,Gf1,P,'sd',2); 
% title('Zoom in Nyquist Diagram') ; 
% legend('Parametric','Non parametric')
% axis equal
% exportgraphics(gcf,'Figures/NyquistData1.png','Resolution',300)

%%

load('data2.mat')
G2 = bj(data2, [5, 5, 5, 5, 1]) ; 
Gf2 = spa(data2, 400) ; 


%% Plot Nyquist 

% P = nyquistoptions; 
% P.ConfidenceRegionDisplaySpacing=1; 
% P.ShowFullContour='off';
% figure();
% subplot(121) 
% nyquistplot(G2,Gf2,P,'sd',2); 
% subplot(122)
% P.Xlim=[45, 50]; P.Ylim=[-30, 5]; 
% h=nyquistplot(G2,Gf2,P,'sd',2); 
% title('Zoom in Nyquist Diagram') ; 
% axis equal

%%

load('data3.mat')
G3 = bj(data3, [5, 5, 5, 5, 1]) ; 
Gf3 = spa(data3, 400) ; 

%% Plot Nyquist 
% P = nyquistoptions; 
% P.ConfidenceRegionDisplaySpacing=1; 
% P.ShowFullContour='off';
% figure();
% subplot(121) 
% nyquistplot(G3,Gf3,P,'sd',2); 
% subplot(122)
% P.Xlim=[45, 50]; P.Ylim=[-30, 5]; 
% h=nyquistplot(G3,Gf3,P,'sd',2); 
% title('Zoom in Nyquist Diagram') ; 
% axis equal

%%

load('data4.mat')
G4 = bj(data4, [5, 5, 5, 5, 1]) ; 
Gf4 = spa(data4, 400) ; 

%% Plot Nyquist 
% P = nyquistoptions; 
% P.ConfidenceRegionDisplaySpacing=1; 
% P.ShowFullContour='off';
% figure();
% subplot(121) 
% nyquistplot(G4,Gf4,P,'sd',2); 
% subplot(122)
% P.Xlim=[45, 50]; P.Ylim=[-30, 5]; 
% h=nyquistplot(G4,Gf4,P,'sd',2); 
% title('Zoom in Nyquist Diagram') ; 
% axis equal


%% compute the filters 
Gmm = stack(1, G1, G2, G3, G4) ; 

Gnom1 = G1 ; 
[Gu1, info1] = ucover(Gmm, Gnom1, 1) ;
W2_1=info1.W1opt ; 

Gnom2 = G2 ; 
[Gu2, info2] = ucover(Gmm, Gnom2, 1) ;
W2_2=info2.W1opt ; 

Gnom3 = G3 ; 
[Gu3, info3] = ucover(Gmm, Gnom3, 1) ;
W2_3=info3.W1opt ; 

Gnom4 = G4 ; 
[Gu4, info4] = ucover(Gmm, Gnom4, 1) ;
W2_4=info4.W1opt ; 

%% plot the magnitude of the filters
% figure()
% bodemag(W2_1)
% hold on ; grid on 
% bodemag(W2_2)
% bodemag(W2_3)
% bodemag(W2_4)
% legend('W2_1','W2_2','W2_3','W2_4') 
%% compute the norm of the filters 
n1 = norm(W2_1,Inf)
n2 = norm(W2_2,Inf)
n3 = norm(W2_3,Inf)
n4 = norm(W2_4,Inf)

%% relative error with filters 

figure()
err1 = Gmm/Gnom1 -1 ; 
bodemag(err1,W2_1) ; grid on ; legend('Gmm/Gnom1 - 1','W_2,1')
% exportgraphics(gcf,'Figures/Err1.png','Resolution',300)

figure()
err2 = Gmm/Gnom2 -1 ; 
bodemag(err2,W2_2) ; grid on; legend('Gmm/Gnom2 - 1','W_2,2')
% exportgraphics(gcf,'Figures/Err2.png','Resolution',300)

figure()
err3 = Gmm/Gnom3 -1 ; 
bodemag(err3,W2_3) ; grid on; legend('Gmm/Gnom3 - 1','W_2,3')
% exportgraphics(gcf,'Figures/Err3.png','Resolution',300)

figure()
err4 = Gmm/Gnom4 -1 ; 
bodemag(err4,W2_4) ; grid on; legend('Gmm/Gnom4 - 1','W_2,4')
% exportgraphics(gcf,'Figures/Err4.png','Resolution',300)


clc

