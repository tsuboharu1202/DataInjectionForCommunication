function sol = G11_merge(n,m,T)
% G10 は 双対Λの対称性 Λ1 - Λ1' = 0
n1_2 = (3*n+m)*(3*n+m);
sol = struct();
sol.dL = sparse(n1_2,n*m);
sol.dY = sparse(n1_2,n*n);
sol.dAlpha = sparse(n1_2,1);
sol.dBeta = sparse(n1_2,1);
sol.dtDelta = sparse(n1_2,1);


sol.dLambda1 = speye(n1_2) - implicit.helper.commutation(3*n+m,3*n+m);
sol.dLambda3 = sparse(n1_2,(n+m)*(n+m));
sol.Lambda_Alpha = sparse(n1_2,1);
sol.Lambda_Beta = sparse(n1_2,1);
sol.Lambda_tDelta = sparse(n1_2,1);
sol.Lambda_Y = sparse(n1_2,n*n);
sol.Data = sparse(n1_2,T*(2*n+m));

sol.G11_row_without_Data = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta, sol.Lambda_Y];
end