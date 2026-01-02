function sol = G1_merge(n,m, B,T)
% G1 は Lagrangian をLで変微分したもの
nm = n*m;
sol = struct();
sol.dL = sparse(nm,nm);
sol.dY = sparse(nm,n*n);
sol.dAlpha = sparse(nm,1);
sol.dBeta = sparse(nm,1);
sol.dtDelta = sparse(nm,1);

sol.dLambda1 = (-implicit.helper.dF1_dL(n,m, B))';
sol.dLambda3 = (-implicit.helper.dF3_dL(n,m))';
sol.Lambda_Alpha = sparse(nm,1);
sol.Lambda_Beta = sparse(nm,1);
sol.Lambda_tDelta = sparse(nm,1);
sol.Lambda_Y = sparse(nm,n*n);
sol.Data = sparse(nm,T*(2*n+m));

sol.G1_row_without_Data = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta, sol.Lambda_Y];
end