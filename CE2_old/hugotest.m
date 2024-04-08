clear ; close ; clc ;

load('data1.mat')

G1 = bj(data1, [5, 5, 5, 5, 1]) ;  %parametric
Gf1 = spa(data1, 400) ;             %non parametric

%% Plot Nyquist 
% 
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

compare( G1 , data1 )

%%

load('data1.mat', 'tt8');

%%

load data1.mat tt8
tt8e = tt8(1:250,:);
tt8v = tt8(251:end,:);

