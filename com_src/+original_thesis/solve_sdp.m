function [sol, K, Y, L, diagnostics] = solve_sdp(data, opts)


if nargin < 2 || isempty(opts), opts = struct(); end
if ~isfield(opts,'verbose'),   opts.verbose   = 0;     end
if ~isfield(opts,'solver'),    opts.solver    = 'mosek'; end

% -------------------------
% Unpack data
% -------------------------
B = data.B;
Z     = data.Z;     % n×T
Xm     = data.X;    % n×T
Um     = data.U;    % m×T  (m=1 でも可)


Phi11  = data.Phi11;
Phi12  = data.Phi12;
Phi22  = data.Phi22;

[n, T] = size(Z);
m = size(Um,1);

% -------------------------
% Decision variables
% -------------------------
% 論文でのXはここではXmと混同しないためにLとする
Y      = sdpvar(n,n,'symmetric');
L      = sdpvar(m,n,'full');
alpha  = sdpvar(1,1);
beta   = sdpvar(1,1);
tDelta = sdpvar(1,1);     % = δ²

% -------------------------
G = [ eye(n)  ,  Z-B*Um ;           % 1: n
    zeros(n,n) , -Xm ;           % 2: n
    zeros(n,n) , zeros(n,T) ;           % 3: m
    zeros(m,n+T)];                 % 7: m

Phi  = [Phi11 Phi12; Phi12' Phi22];

% -------------------------
% LMI blocks (25)(26)(27)
% -------------------------
[F1, F2, F3] = original_thesis.build_lmi_blocks(Y,L,alpha,beta,tDelta,G,Phi,data);

% -------------------------
% Constraints
% -------------------------

tolerance = 1e-8;

const_1 = [F1 - F2 >= 0];
const_2 = [F3 >= tolerance*eye(n+m)];


constr  = [];
constr  = [constr, const_1];
constr  = [constr, const_2];
constr  = [constr, alpha >= 0, beta >= tolerance, tDelta >= tolerance, Y >= tolerance*eye(n), (1 - tolerance)^2 >= tDelta];

% Objective: maximize δ  <=>  minimize tDelta
obj = -tDelta;

params = sdpsettings('solver', opts.solver, 'verbose', opts.verbose);
diagnostics = optimize(constr, obj, params);

if diagnostics.problem ~= 0
    warning('Optimization problem status: %d - %s', diagnostics.problem, diagnostics.info);
end

% -------------------------
% Output pack
% -------------------------
sol.status    = diagnostics.problem;
sol.Y         = value(Y);
sol.L         = value(L);
sol.alpha     = value(alpha);
sol.beta      = value(beta);
sol.tDelta    = value(tDelta);
% tDelta = δ^{2}なので、δ = (tDelta)^(1/2)
delta = sqrt(value(tDelta));
sol.rho = (1-delta)/(1+delta);  % Theorem 3: ρ = (1-δ)/(1+δ)
sol.objective = value(obj);

K = sol.L/sol.Y;
sol.K = K;
end
