function Kdd = toTF(controller)
%toTF Convert controller structure to Transfer Function
%   
%   ------------------------------------------------------------------------
%   Copyright 2021 Philippe Schuchert, DDMAC, EPFL
% 

Kdd = tf( ...
    conv(controller.num, controller.Fx), ...
    conv(controller.den, controller.Fy), ...
    controller.Ts);
end