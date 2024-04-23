function [controller, sol] = datadriven(system, obj, cons, params, verbose, solver)
%DATADRIVEN Synthesis of datadriven controller  
%   This function iteratively solve the datadriven controller synthesis problem.
% 
%   Outputs are ...
% 
%   ------------------------------------------------------------------------
%   Copyright 2024 Vaibhav Gupta, DDMAC, EPFL (MIT License)
%

arguments
    system
    obj
    cons
    params
    verbose (1, 1) logical = true
    solver  (1, 1) string = ""
end

%% Settings and Auxilary functions

% Add solver to the parameter list
params.solver = solver;

% Verbosity Printer
printer.div = @(varargin) fprintf("  |------|------------|------------|-------------|------------|\n");
printer.title = @(varargin) fprintf("  | Iter |   Slack    |  Objective |   Decrease  | Solve Time |\n");
printer.iter = @(varargin) fprintf("  |  %03d | %10.4e | %10.4e | %+11.4e | %10.2f |\n", varargin{:});
printer.head = @(varargin) cellfun(@(x) x(varargin), {printer.div, printer.title, printer.div});

% Disable printer if verbosity is false
if ~verbose
    printer = structfun(@(x) @(varargin) fprintf(""), printer, "UniformOutput",false);
end

% `iff` function gives 
%   * valT , if condition is true 
%   * valF , otherwise
iff = @(condition, valT, valF) subsref({valF, valT}, struct("type", "{}", "subs", {{condition+1}})); 

%% Iterations
fprintf("INFO: Solving for optimal controller in datadriven framework...\n");
printer.head(); % Print table head

prev_obj = [];
sol = [];

solveTime = 0;
for iter = 1:params.maxIter
    switch params.solver
        case 'fusion'
            [tmp_controller, tmp_sol, diagnostics] = solveddMOSEK(system,obj,cons,params,sol);
        otherwise
            [tmp_controller, tmp_sol, diagnostics] = solveddYALMIP(system,obj,cons,params,sol);     
    end
    solveTime = solveTime + diagnostics.solvertime; % Does not include yalmiptime!

    printer.iter( ...
        iter, ...                                                                   # Iteration
        iff(~tmp_sol.satisfyConstraints, tmp_sol.slack, []), ...                    # Slack
        iff( tmp_sol.satisfyConstraints, tmp_sol.obj, []), ...                      # Objective
        iff(tmp_sol.satisfyConstraints, tmp_sol.obj, tmp_sol.slack) - prev_obj, ... # Decrease in objective/slack
        solveTime ...                                                               # Running Solve Time
        );

    % If objective increases, QUIT!!!
    if tmp_sol.obj - prev_obj > 0
        break
    end

    % Update sol and controller
    sol = tmp_sol;
    system.controller = tmp_controller;
    
    % Converged solution check
    if sol.satisfyConstraints && any(sol.obj-prev_obj >= -params.tol)
        break
    end

    % Update prev. objective
    prev_obj = sol.obj;

    % Constraint satisfaction check
    if ~sol.satisfyConstraints && sol.slack < 1e-4
        sol.satisfyConstraints  = 1;
        prev_obj = [];
    end
end

printer.div();  % Print ending divider

%% Returns
controller = toTF(system.controller);
sol = rmfield(sol,'satisfyConstraints');

end