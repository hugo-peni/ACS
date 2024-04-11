% Load empty structure
[SYS, OBJ, CON, PAR] = datadriven.utils.emptyStruct();

% ------------------------------------------------------------------------------
%   Initial controller
% ------------------------------------------------------------------------------
Kc = ...                        % Initial controller
[num, den] = tfdata(Kc, 'v');   % Extract numerator and denominator

den(order + 1) = 0; % Zero padding to have same order as desired controller
num(order + 1) = 0; % Zero padding to have same order as desired controller

% Fixed parts of the controller
%   NOTE: Initial controller should contain the fixed parts too!
Fy = 1;                 % Fixed part of denominator as polynomial
den = deconv(den, Fy);  % Remove fixed part of denominator

Fx = 1;                 % Fixed part of numerator as polynomial
num = deconv(num, Fx);  % Remove fixed part of numerator

SYS.controller.num = num;
SYS.controller.den = den;
SYS.controller.Ts = Ts;
SYS.controller.Fx = Fx;
SYS.controller.Fy = Fy;

% ------------------------------------------------------------------------------
%   Nominal system(s)
% ------------------------------------------------------------------------------
% Systems should be LTI systems (`ss`, `tf`, `frd`, ...)
SYS.model = ...

% ------------------------------------------------------------------------------
%   Frequencies for controller synthesis
% ------------------------------------------------------------------------------
SYS.W = ...

% ------------------------------------------------------------------------------
%   Filters for objectives
% ------------------------------------------------------------------------------
% Filter should be LTI systems (`ss`, `tf`, `frd`, ...)
% For unused objectives, set filters to []
OBJ.oinf.W1 = [];   % ║W1 S║∞
OBJ.oinf.W2 = [];   % ║W2 T║∞
OBJ.oinf.W3 = [];   % ║W3 U║∞
OBJ.oinf.W4 = [];   % ║W4 V║∞

OBJ.o2.W1 = [];     % ║W1 S║₂
OBJ.o2.W2 = [];     % ║W2 T║₂
OBJ.o2.W3 = [];     % ║W3 U║₂
OBJ.o2.W4 = [];     % ║W4 V║₂

OBJ.LSinf.Ld = [];  % ║W (Ld - G K)║∞
OBJ.LSinf.W  = [];

OBJ.LS2.Ld = [];    % ║W (Ld - G K)║₂
OBJ.LS2.W  = [];

% ------------------------------------------------------------------------------
%   Filters for constraints
% ------------------------------------------------------------------------------
% Filter should be LTI systems (`ss`, `tf`, `frd`, ...)
% For unused constraints, set filters to []
CON.W1 = [];    % ║W1 S║∞ ≤ 1
CON.W2 = [];    % ║W2 T║∞ ≤ 1
CON.W3 = [];    % ║W3 U║∞ ≤ 1
CON.W4 = [];    % ║W4 V║∞ ≤ 1

% ------------------------------------------------------------------------------
%   Optimisation parameters
% ------------------------------------------------------------------------------
PAR.tol = 1e-4;     % Numerical tolerance for convergence
PAR.maxIter = 100;  % Maximum number of allowed iterations

verbosity = true;   % To print controller synthesis iterations
solver = "";        % Solver to use for optimisation ("mosek", "sedumi", ...)

% ------------------------------------------------------------------------------
%   Solve the datadriven controller synthesis problem
% ------------------------------------------------------------------------------
[K, sol] = datadriven.datadriven(SYS, OBJ, CON, PAR, verbosity, solver);

