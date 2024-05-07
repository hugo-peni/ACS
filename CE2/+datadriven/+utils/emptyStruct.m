function [SYS, OBJ, CON, PAR] = emptyStruct()
%emptyStruct Creates the required empty structures for datadriven
% 
%   ------------------------------------------------------------------------
%   Copyright 2021 Philippe Schuchert, DDMAC, EPFL
%   Copyright 2024 Vaibhav Gupta, DDMAC, EPFL (MIT License)
% 

% Controller
CTRL = struct('num',[], 'den',[], 'Ts',[], 'Fx',[], 'Fy',[]);

% System and controller
SYS = struct('model',[], 'W',[], 'controller', CTRL);

% Parameters
PAR = struct( ...
    'tol', 1e-6, ...
    'maxIter', 100, ...
    'radius', 1, ...
    'robustNyquist', true, ...
    'robustEllipsoid', []);

% Objectives
% 
%   ║W1 S║
%   ║W2 T║ 
%   ║W3 U║
%   ║W4 V║₂/∞
% 
%   ║W (Ld - G K)║₂/∞
OBJ  = struct( ...
    'oinf' , struct('W1',[], 'W2',[], 'W3',[], 'W4',[]), ...
    'o2'   , struct('W1',[], 'W2',[], 'W3',[], 'W4',[]), ...
    'LSinf', struct('W' ,[], 'Ld',[]), ...
    'LS2'  , struct('W' ,[], 'Ld',[]));

% Constraint Defination
% 
%   ║W1 S║ ≤ 1
%   ║W2 T║ 
%   ║W3 U║
%   ║W4 V║∞
% 
CON = struct('W1',[], 'W2',[], 'W3',[], 'W4',[]);

end