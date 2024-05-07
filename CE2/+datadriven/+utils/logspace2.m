function w = logspace2(a,b,c)
% 
%   ------------------------------------------------------------------------
%   Copyright 2021 Philippe Schuchert, DDMAC, EPFL
% 
w = logspace(log10(a),log10(b),c);
w(end) = b;
end
