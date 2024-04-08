function parameters = datadrivenOptions(varargin)
% datadrivenOptions(name1,value1, name2, value2, ...)
% Structure of parameters used in datadriven command.
% --------------------------------------------------------------------
%     name          ->          description
% --------------------------------------------------------------------
% 'K'               -> intitial stabilizing controller
% 'order'           -> order of the final controller
% 'force_integrator'-> forces an integrator in the final controller
% 'W'               -> Frequency grid where the problem is solved
% 'max_iteration'   -> number of iteration before stopping
% 'sdpsettings'     -> sdpsettings structure for YALMIP
% --------------------------------------------------------------------
%
% Example use:
% options = datadrivenOptions('Kc',Kinit, 'order', orderK, ...
%                             'W',W, ...
%                             'force_integrator',true, ...
%                             'max_iteration',20);

parameters = struct('Kc',[], 'force_integrator',false, 'W',[], 'max_iteration',10, 'sdpsettings',sdpsettings('verbose',0));

tokens = {};
if mod(nargin,2) ~= 0
    error('arguments must come in pairs');
end

for ii = 1 : 2: nargin
    key  = varargin{ii};
    value= varargin{ii+1};
    switch lower(deblank(key))
        case 'order'
            tokens{end+1} = struct('key', key, 'value', value);
        otherwise
            if isfield(parameters, key)
                parameters.(key) = value;
            else
                warning('Parameters %s not recoginosed', key)
            end
    end
end
for token = [tokens{:}]
    % only process theses values after the rest has been set
    switch lower(token.key)
        case 'order'
            try
                [a,~,Ts] = tfdata(parameters.Kc, 'v');
            catch
                error('Controller should be a state-space, zpk, or transfer function with correct sampling time')
            end
            if length(a)-1 < token.value
                order_ =  token.value -(length(a)-1);
                A = triu(-0.5*ones(order_)) + triu(-0.5*ones(order_),1); % Padding for correct order
                B = ones(order_,1);
                C = zeros(1,order_);
                parameters.Kc = balreal(parameters.Kc + 0*ss(A,B,C,0,Ts)); % augment order. 
            end
            if length(a)-1 > token.value
                error('Initial controller can not be of higher order than final controller')
            end
    end
end
end