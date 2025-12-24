function sol = G9_merge(n,m,T,beta,Lambda_beta)
% G9 は 相補性条件における βΛ_beta = 0
sol = struct();
sol.dL = sparse(1,n*m);
sol.dY = sparse(1,n*n);
sol.dAlpha = sparse(1,1);
sol.dBeta = Lambda_beta;
sol.dtDelta = sparse(1,1);


sol.dLambda1 = sparse(1,(3*n+m)*(3*n+m));
sol.dLambda3 = sparse(1,(n+m)*(n+m));
sol.Lambda_Alpha = sparse(1,1);
sol.Lambda_Beta = beta;
sol.Lambda_tDelta = sparse(1,1);
sol.Data = sparse(1,T*(2*n+m));

sol.G9_row_without_Data = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta];
end