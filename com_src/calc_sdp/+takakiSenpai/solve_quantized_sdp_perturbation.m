% function [sol, K, Y, L, diagnostics] = solve_quantized_sdp_perturbation(data, opts)
%     % Solve SDP (25)(26)(27) and maximize delta via tDelta = δ^{-2} minimization.
    

%     if nargin < 2 || isempty(opts), opts = struct(); end
%     if ~isfield(opts,'verbose'),   opts.verbose   = 0;     end
%     if ~isfield(opts,'solver'),    opts.solver    = 'sedumi'; end

%     % -------------------------
%     % Unpack data
%     % -------------------------
%     Z     = data.Z;     % n×T
%     Xm     = data.X;    % n×T
%     Um     = data.U;    % m×T  (m=1 でも可)

%     H      = data.H;
%     E1     = data.E1;
%     E2     = data.E2;
%     n2     = size(E2,1);      

%     Phi11  = data.Phi11;
%     Phi12  = data.Phi12;
%     Phi22  = data.Phi22;

%     [n, T] = size(Z);
%     m = size(Um,1);

%     % -------------------------
%     % Decision variables
%     % -------------------------
%     Y      = sdpvar(n,n,'symmetric');
%     L      = sdpvar(m,n,'full');
%     alpha  = sdpvar(1,1);
%     beta   = sdpvar(1,1);
%     tTau   = sdpvar(1,1);     % = τ^{-2}
%     tDelta = sdpvar(1,1);     % = δ^{-2}

%     % -------------------------
%     G = [ eye(n)  ,  Z ;           % 1: n
%       zeros(n,n) , -Xm ;           % 2: n
%       zeros(m,n) , -Um ;           % 3: m
%       zeros(m,n+T) ;                  % 4: m
%       zeros(n,n+T) ;                  % 5: n
%       zeros(n2,n+T);                  % 6: n2
%       zeros(m,n+T) ];                 % 7: m

%     Phi  = [Phi11 Phi12; Phi12' Phi22];

%     % -------------------------
%     % LMI blocks (25)(26)(27)
%     % -------------------------
%     [M25, M26, M27] = build_lmi_blocks_perturbation(Y, L, alpha, beta, tTau, tDelta, H, E1, E2, Phi,Z,Xm,Um);

%     % -------------------------
%     % Constraints
%     % -------------------------
%     if ~isfield(opts,'psd_slack'), opts.psd_slack = 1e-10; end
    
%     constr  = [];
%     constr  = [constr, Y >= opts.psd_slack*eye(n)];
%     constr  = [constr, alpha >= 0, beta >= opts.psd_slack, tTau >= opts.psd_slack, tDelta >= opts.psd_slack + 1];
%     constr  = [constr, M25 >= 0];
    
%     constr = [constr, M26 - opts.psd_slack*eye(size(M26,1)) >= 0];
%     constr = [constr, M27 - opts.psd_slack*eye(size(M27,1)) >= 0];

%     % Objective: maximize δ  <=>  minimize tDelta
%     obj = tDelta;

%     % -------------------------
%     % Solve
%     % -------------------------
%     params        = sdpsettings('solver', opts.solver, 'verbose', opts.verbose);
%     diagnostics = optimize(constr, obj, params);
%     fprintf('status = %d (%s)\n', diagnostics.problem, yalmiperror(diagnostics.problem));
    
%     % if diagnostics.problem ~= 0
%     %     warning('SDP failed: %s', yalmiperror(diagnostics.problem));
%     %     sol.status = diagnostics.problem;
%     %     sol.tDelta = NaN;
%     %     K = NaN; Y = NaN; L = NaN;
%     %     return;
%     % end

%     % -------------------------
%     % Output pack
%     % -------------------------
%     sol.status    = diagnostics.problem;
%     sol.Y         = value(Y);
%     sol.L         = value(L);
%     sol.alpha     = value(alpha);
%     sol.beta      = value(beta);
%     sol.tTau      = value(tTau);
%     sol.tDelta    = value(tDelta);
%     delta = (value(tDelta))^(-1/2);
%     sol.rho = (1+delta)/(1-delta);
%     sol.objective = value(obj);

%     try
%         K = recover_gain(sol.Y, sol.L);
%     catch
%         K = NaN;
%     end
% end
