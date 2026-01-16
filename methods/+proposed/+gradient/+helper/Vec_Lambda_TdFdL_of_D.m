function v = Vec_Lambda_TdFdL_of_D(n,m,T,D,Lambda)
% v = vec( grad_L <Lambda, F(L)> ) as a function of D=[Z;X;U]
% Here grad w.r.t L is:
%   G(D) = X'*(Lambda11+Lambda22) + 2*Z'*Lambda21 + 2*U'*Lambda41
% where blocks correspond to F1 structure (2n+2m).

Z = D(1:n, :);
X = D(n+1:2*n, :);
U = D(2*n+1:2*n+m, :);

Lambda = 0.5*(Lambda + Lambda'); % safety

Lambda11 = Lambda(1:n, 1:n);
Lambda21 = Lambda(n+1:2*n, 1:n);
Lambda22 = Lambda(n+1:2*n, n+1:2*n);
Lambda41 = Lambda(2*n+m+1:2*n+2*m, 1:n);

G = X'*(Lambda11 + Lambda22) + 2*Z'*Lambda21 + 2*U'*Lambda41; % (T x n)
v = G(:);
end
