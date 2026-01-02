function sol = G2_merge(n,m, B,T)
% G2 は Lagrangian をYで変微分したもの
nn = n*n;
nm = n*m;
sol = struct();
sol.dL = sparse(nn,nm);
sol.dY = sparse(nn,nn);
sol.dAlpha = sparse(nn,1);
sol.dBeta = sparse(nn,1);
sol.dtDelta = sparse(nn,1);

sol.dLambda1 = (-implicit.helper.dF1_dY(n,m))';
sol.dLambda3 = (-implicit.helper.dF3_dY(n,m))';
sol.Lambda_Alpha = sparse(nn,1);
sol.Lambda_Beta = sparse(nn,1);
sol.Lambda_tDelta = sparse(nn,1);
sol.Lambda_Y = speye(nn);
sol.Data = sparse(nn,T*(2*n+m));

sol.G2_row_without_Data  = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta, sol.Lambda_Y];
end