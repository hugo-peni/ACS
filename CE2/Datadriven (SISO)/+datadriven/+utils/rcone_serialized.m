function rcone = rcone_serialized(x,y,z)
%RCONE_SERIALIZED Seralised version of rotated cone
%   Cone command in YALMIP is not serialised. This function provides the
%   serialisation.
% 
%   ------------------------------------------------------------------------
%   Copyright 2021 Philippe Schuchert, DDMAC, EPFL
%   Copyright 2024 Vaibhav Gupta, DDMAC, EPFL (MIT License)
% 

if ~isempty(z)
    % ‖z‖^2 < 2 x y
    % x + y > 0
    v = 1/sqrt(2);
    T = blkdiag([v,v; v,-v], eye(2*size(z,2)));
    Z = T\[x.';y.';real(z).';imag(z).'];
    rcone = cone(Z);
else
    rcone = [];
end
end