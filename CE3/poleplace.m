function [Rp,Sp]=poleplace(B,A,Hr,Hs,P) 

    %Compute the convolution and orders
    A_p = conv(A,Hs) ; B_p = conv(B,Hr) ; 
    nA = length(A_p) - 1 ; nB = length(B_p) - 1 ; 
    
    d = find(B_p==0, 1, 'last' ) - 1   ;
    nS = nB + d - 1 ;
    
    %Compute Sylverster Matrix
    M = zeros(nA + nB + d, nA + nB + d) ; 
    for ind=1:nB+d
        M(ind:nA+ind , ind) = A_p;
    end
    for ind=1:nA
        M(ind+d:ind+nB+d , ind+nB+d) = B_p ;
    end
    
    %Solve equation and retrieve R and S
    p=zeros(nA+nB+d,1) ; Ps = length(P) ;
    p(1:Ps) = transp(P) ;
    
    x = M \ p ;
    S = x(1:nS+1)' ;
    R = x(nS+2:end)' ; 

    Sp = conv(Hs,S) ; 
    Rp = conv(Hr,R) ; 


end

