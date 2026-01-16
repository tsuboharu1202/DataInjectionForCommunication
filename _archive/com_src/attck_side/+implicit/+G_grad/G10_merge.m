function sol = G10_merge(n,m,T,tDelta,Lambda_tDelta)
% G10 は 相補性条件における tΔΛ_tΔ = 0
sol = struct();
sol.dL = sparse(1,n*m);
sol.dY = sparse(1,n*n);
sol.dAlpha = sparse(1,1);
sol.dBeta = sparse(1,1);
sol.dtDelta = Lambda_tDelta;


sol.dLambda1 = sparse(1,(3*n+m)*(3*n+m));
sol.dLambda3 = sparse(1,(n+m)*(n+m));
sol.Lambda_Alpha = sparse(1,1);
sol.Lambda_Beta = sparse(1,1);
sol.Lambda_tDelta = tDelta;
sol.Lambda_Y = sparse(1,n*n);
sol.Data = sparse(1,T*(2*n+m));

sol.G10_row_without_Data = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta,sol.Lambda_Y];
end