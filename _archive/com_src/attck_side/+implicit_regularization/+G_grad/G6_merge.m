function sol = G6_merge(n,m,T)
% G6 は Lambda_P = Lambda_P'
n1 = 2*n+2*m;
sol = struct();
sol.dL = sparse(n*n,n*T);
sol.dLambda = sparse(n*n,n1*n1);
sol.dtDelta = sparse(n*n,1);
sol.dLambda_P = speye(n*n)-implicit.helper.commutation(n,n);

sol.Data = sparse(n*n,T*(2*n+m));
sol.dLambda3 = sparse(n*n, n*n);
sol.G6_row_without_Data = [sol.dtDelta, sol.dL, sol.dLambda,sol.dLambda_P, sol.dLambda3];
end