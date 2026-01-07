function test_dVec_Lambda_TdFdL_dD_FD1(n,m,T,D,Lambda)
% FD1本で dVec_Lambda_TdFdL_dD の正しさを確定
% D: (2n+m) x T

fprintf('\n=== dVec(Lambda^T dF/dL)/dD FD1 check (single direction) ===\n');

J = implicit_regularization.helper.dVec_Lambda_TdFdL_dD(n,m,T,Lambda); % (T*n) x (T*(2n+m))

v0 = implicit_regularization.helper.Vec_Lambda_TdFdL_of_D(n,m,T,D,Lambda);

E = randn(size(D));
e = E(:);

eps_list = [1e-6, 1e-7, 1e-8];

for k = 1:numel(eps_list)
    eps = eps_list(k);
    
    v1 = implicit_regularization.helper.Vec_Lambda_TdFdL_of_D(n,m,T, D + eps*E, Lambda);
    
    lhs = (v1 - v0)/eps; % FD
    rhs = J * e;         % analytic
    
    relerr = norm(lhs - rhs) / max(1, norm(lhs));
    fprintf('  eps=%.1e: relerr=%.3e\n', eps, relerr);
end

fprintf('=== done ===\n');
end
