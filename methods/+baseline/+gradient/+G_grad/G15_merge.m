function sol = G15_merge(n,m,T,Lambda_Y)
% G10 は 双対Λの対称性 Lambda_Y - Lambda_Y' = 0
n1_2 = n*n;
sol = struct();
sol.dL = sparse(n1_2,n*m);
sol.dY = sparse(n1_2,n*n);
sol.dAlpha = sparse(n1_2,1);
sol.dBeta = sparse(n1_2,1);
sol.dtDelta = sparse(n1_2,1);


sol.dLambda1 = sparse(n1_2,(3*n+m)*(3*n+m));
sol.dLambda3 = sparse(n1_2,(n+m)*(n+m));
sol.Lambda_Alpha = sparse(n1_2,1);
sol.Lambda_Beta = sparse(n1_2,1);
sol.Lambda_tDelta = sparse(n1_2,1);
sol.Lambda_Y = speye(n1_2)-core.helper.commutation(n,n);
sol.Data = sparse(n1_2,T*(2*n+m));

sol.G15_row_without_Data = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta, sol.Lambda_Y];
end