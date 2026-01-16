function test_dPi_dD_FD1(n,m,T,D)
% FD1本で dPi_dD の正しさを確定
% D: (2n+m) x T, ordered as [Z; X; U]

fprintf('\n=== dPi_dD FD1 check (single direction) ===\n');

% analytic Jacobian
Gamma = [D(2*n+1:2*n+m,:); D(n+1:2*n,:)]; % [U;X]
J = helper.dPi_dD(n,m,T,Gamma); % (T^2) x (T*(2n+m))

% base Pi
Pi0 = helper.Pi_of_D(n,m,T,D);
v0  = Pi0(:);

% random direction (single)
E = randn(size(D));
e = E(:);

eps_list = [1e-6, 1e-7, 1e-8];

for k = 1:numel(eps_list)
    eps = eps_list(k);
    
    Pi1 = helper.Pi_of_D(n,m,T, D + eps*E);
    v1  = Pi1(:);
    
    lhs = (v1 - v0)/eps;      % FD
    rhs = J * e;              % analytic
    
    relerr = norm(lhs - rhs) / max(1, norm(lhs));
    fprintf('  eps=%.1e: relerr=%.3e\n', eps, relerr);
end

fprintf('=== done ===\n');
end
