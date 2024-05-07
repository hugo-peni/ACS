function [X,Y] = toXY(controller)
%toXY Convert controller structure to numerator and denominator form
%   K = X/Y
% 
%   ------------------------------------------------------------------------
%   Copyright 2021 Philippe Schuchert, DDMAC, EPFL
% 
X = tf(conv(controller.num,controller.Fx), 1, controller.Ts);
Y = tf(conv(controller.den,controller.Fy), 1, controller.Ts);
end