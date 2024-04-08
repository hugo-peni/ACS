function K_optimal = datadriven_ACS(G, OBJ, CON, PAR)
% Author: Philippe Schuchert
%
% datadriven :  Data-driven mixed-sensitivity controller design. SISO only.
% Computes the optimal controller for the mixed synthesis problem.
% For the feedback loop
%
%                     e       u
%            r --->O--->[ K ]-->[ G ]---+---> y
%                - |                    |
%                  +<-------------------+
%
%    the mixed-sensitivity design seeks a controller K that minimizes
%
%                   || OBJ.W1*S ||
%             min   || OBJ.W2*U ||
%                   || OBJ.W3*T ||_norm
%        subject to
%        ||  CON.W1*S  ||_\infty ≤ 1
%        ||  CON.W2*U  ||_\infty ≤ 1
%        ||  CON.W3*T  ||_\infty ≤ 1
%
%    where
%      * S = 1/(I+G*K) is the sensitivity function
%      * U is the transfer from r to u (control effort)
%      * T = I-S = G*K/(I+G*K) is the complementary sensitivity
%      * norm is the norm specified in options (either 'two' or 'inf')
% If multi-model uncertainity is used, use G = cat(1,G1,G2,..);



% K = X/Y, where  
% [X;Y] = C/(zI - A)*B + D
% and C, D are optimization parameters. 
% Parameters A, B are recomputed every iteration 

C = sdpvar(2, order(PAR.Kc), 'full');
D = sdpvar(2, 1, 'full');

if ~isa(PAR,'struct')
    error('Last argument should be a datadrivenOptions structure')
end

if any(isnan(freqresp(G, PAR.W)),'all')
    error('NaN when evaluating the frequency response of G. Check if frequencies appropriately chosen')
end


fprintf('\n\n -------------------------------------------------------------------\n')
fprintf(' Starting new optimization \n')


% Does the inital controller satisfy all H infinity constraints?
fprintf(' -------------------------------------------------------------------\n')
fprintf(' Searching for initial controller satisfying all Hinf constraints\n')
fprintf(' -------------------------------------------------------------------\n')

PAR = find_initial_controller(G, OBJ, CON, PAR, C, D);

fprintf(' -------------------------------------------------------------------\n')
fprintf(' Optimizing over objective\n')
fprintf(' -------------------------------------------------------------------\n')

PAR = find_optimal_controller(G, OBJ, CON, PAR, C, D);

K_optimal = PAR.Kc;
end

%%  Main functions

function PAR = find_initial_controller(G, ~, CON, PAR, C, D)
W = PAR.W;
upper_bound = sdpvar(sdpvar);
for iter = 1 : PAR.max_iteration
    fprintf(' Iteration %02d',iter)
    
    % Frequency response controller
    [X,Y,Xc,Yc,A,B] = controllers_data(PAR, W, D, C);
    
    if PAR.force_integrator
        cosntraints = C(1,:)/(eye(length(A))-A)*B + D(1) == 0;
    else
        cosntraints = [];
    end
    
    for mod = 1 : length(G)
        Gf = resp(G(:,:,mod), W);
        
        P = Y + Gf.*X;
        Pc = Yc + Gf.*Xc;
        PHI = conj(P).*Pc + conj(Pc).*P - conj(Pc).*Pc;
        
        % Additional stability constraint, for safety.
        cosntraints = [cosntraints,  polygonalStability(P,Pc)];
        
        % Soft constraints
        GAM = 1 + upper_bound*ones(length(W), 1);
        if ~isempty(CON.W1)
            cosntraints = [cosntraints, rcone_serialized(PHI, GAM, resp(CON.W1,W).*Y) ];
        end
        if ~isempty(CON.W2)
            cosntraints = [cosntraints, rcone_serialized(PHI, GAM, resp(CON.W2,W).*X) ];
        end
        if ~isempty(CON.W3)
            cosntraints = [cosntraints, rcone_serialized(PHI, GAM, resp(CON.W3,W).*Gf.*X) ];
        end
    end
    
    JOB = optimize([upper_bound >=0, cosntraints], upper_bound, PAR.sdpsettings);
    fprintf(' | Slack in constraints: %3.4e | %s\n',sqrt(double(upper_bound)), JOB.info)
    
    if double(upper_bound) < 1e-8
        fprintf('                Slack ≤ 1E-4. Can proceed to objective\n')
        break
    end
    
    YX = ss(A, B, double(C), double(D), PAR.Kc.Ts);
    PAR.Kc =  reform(YX);
    
    if iter == PAR.max_iteration
        error('   Could not find appropirate initial controller. Aborting\n')
    else
    end
end
end

function PAR = find_optimal_controller(G,OBJ,CON,PAR,C,D)
W = PAR.W;
for iter = 1 : PAR.max_iteration
    fprintf(' Iteration %02d',iter)
    
    % Controller data
    [X, Y, Xc, Yc, A, B] = controllers_data(PAR, W, D, C);
    
    if PAR.force_integrator
        constraints = C(1,:)/(eye(length(A))*(PAR.Kc.Ts>0)-A)*B + D(1) == 0;
    else
        constraints = [];
    end
    
    if  strcmpi(OBJ.norm,'two') 
        % H2 optimization
        gamma = sdpvar(length(W), length(G), 'full');
        local_h2_norm = sum(gamma.*kron([diff(W(:));0]+[W(1);diff(W(:))], ones(1, length(G)))); % corresponds  to trapz(W,gamma(:,mod))
        obj = max(local_h2_norm)/(2*pi); % optimize worste case H2 norm
        if PAR.Kc.Ts
            obj = obj*PAR.Kc.Ts;
        end
    elseif strcmpi(OBJ.norm,'inf') 
        % Hinf optimization
        gamma = sdpvar*ones(numel(W(:)),length(G));
        obj = gamma(1);
    else
        error('Wrong name for norm')
    end
    
    for mod = 1 : length(G)
        Gf = resp(G(:,:,mod),W);
        
        P = Y + Gf.*X;
        Pc = Yc + Gf.*Xc;
        PHI = conj(P).*Pc + conj(conj(P).*Pc) - conj(Pc).*Pc;
        
        % Additional stability constraint, for safety.
        constraints = [constraints,  polygonalStability(P,Pc)];
        
        % Objective
        F_ = [];
        if ~isempty(OBJ.W1)
            F_ = resp(OBJ.W1,W).*Y;
        end
        if ~isempty(OBJ.W2)
            F_ = [F_,resp(OBJ.W2,W).*X];
        end
        if ~isempty(OBJ.W3)
            F_ = [F_,resp(OBJ.W3,W).*Gf.*X ];
        end
        constraints = [constraints, rcone_serialized(PHI, gamma(:,mod), F_)];
        
        % Hinf constraints
        GAM = ones(length(W),1);
        if ~isempty(CON.W1)
            constraints = [constraints, rcone_serialized(PHI, GAM, resp(CON.W1,W).*Y)];
        end
        if ~isempty(CON.W2)
            constraints = [constraints, rcone_serialized(PHI, GAM, resp(CON.W2,W).*X)];
        end
        if ~isempty(CON.W3)
            constraints = [constraints, rcone_serialized(PHI, GAM, resp(CON.W3,W).*Gf.*X)];
        end
        
    end
    JOB = optimize(constraints, obj, PAR.sdpsettings);
    fprintf(' | Objective: %3.4e | %s\n',sqrt(double(obj)), JOB.info)
    
    YX = ss(A,B,double(C),double(D),PAR.Kc.Ts);
    PAR.Kc =  reform(YX);
    
end
fprintf(' Maximum iteration reached. \n')
fprintf(' -------------------------------------------------------------------\n')

end

%% Helper functions

function [X,Y,Xc,Yc,A,B] = controllers_data(PAR,W,D,C)
% Get frequency response of controllers (and matrix A, C)
YX = rncf(PAR.Kc);

[A,B,C_,D_] = ssdata(YX);
assign(C,C_);
assign(D,D_);

Z = resp(tf([1,0],[0 1],PAR.Kc.Ts),W);

KER = ssfresp(A,B,[],[],[],Z);
KER = reshape(KER,[length(A), numel(W)]);
YXf = C*KER + D*ones(1,length(W));
Y = YXf(1,:).';
X = YXf(2,:).';
Yc = double(Y);
Xc = double(X);
end

function F = resp(G,W)
% shortcut for freqresp function
if isa(G,'double')
    F = G + 0*W(:);
else
    f = freqresp(G,W);
    F = f(:);
end
end

function rcone = rcone_serialized(x,y,z)
% SISO case [a,b,b',c]>0 is equivalent to a*c > f'*f, a,c > 0.
% This can be implemented efficiently using a second order rotated code
T = blkdiag([1,1;1,-1], eye(2*size(z,2)));
Z = T\[x.';y.';real(z).';imag(z).'];
rcone = cone(Z);
end

function stab = polygonalStability(P,Pc)
% Additional stability constraint, for safety. See "A convex set of robust
% D—stabilizing controllers using Cauchy's argument principle" by Philippe Schuchert
% Here the set D is the unit circle.

r = [conj(Pc(1));Pc;conj(Pc(end))];
PolygonalChain = [conj(P(1));P;conj(P(end))];

for ii = 1 : length(r)-1
    A = [real(r(ii)),imag(r(ii))];
    B = [real(r(ii+1)),imag(r(ii+1))];
    M = B - A;
    t0 = dot(M, -A) / dot(M, M);
    if t0 < 0
        C = A;
    elseif t0 < 1
        C = A + t0 * M;
    else
        C = B;
    end
    n(ii) = C(1) + 1j*C(2);
end
n = n(:)./abs(n(:));

x1_a = real(PolygonalChain(2:end,:).*conj(n));
x1_b = real(PolygonalChain(1:end-1,:).*conj(n));

stab = [x1_b >= 1e-6, x1_a >= 1e-6]; % 1e-6: some tolerance to numerical solver
end

function K = reform(YX)
% Compute minimal realization of K = X/Y; See  "Data-driven fixed-structure
% frequency-based H2 and H∞ controller design" by Philippe Schuchert 
Y = YX(1);
X = YX(2);

K = balreal(ss(X.A-(X.B/Y.D)*Y.C,(X.B/Y.D), X.C - (X.D/Y.D)*Y.C, (X.D/Y.D),YX.Ts));

end