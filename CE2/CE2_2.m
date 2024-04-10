%%

Gnom = G(:,:,6,1) ;
[Gu, info] = ucover(G, Gnom, 4) ;
W2 = info.W1 ;

%%

% Define the weighting filter W1(s) in continuous time
s = tf('s') ; m = 0.505 ; omega_b = 1000 ;
W1s = m * (s+omega_b) / (s+1e-5) ;
Ts = 0.002 ; W1z = c2d(W1s, Ts, 'zoh');


figure(1)
bode( 1 / W1s )
inf_norm = norm( 1 / W1s , inf ) 
inf_norm_db = mag2db( inf_norm )

K = mixsyn( tf(Gnom) , W1z , [] , W2 ) ;

%%

% Define the CL transfer functions

S=feedback(1,G*K) ; T=feedback(G*K,1) ; U=feedback(K,G) ;

