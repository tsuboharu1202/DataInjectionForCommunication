function sol = G3_merge(n,m,T,Lambda1,F2,B,G,Phi)
% G3 は Lagrangian をAlphaで変微分したもの
sol = struct();
sol.dL = sparse(1,n*m);
sol.dY = sparse(1,n*n);
sol.dAlpha = sparse(1,1);
sol.dBeta = sparse(1,1);
sol.dtDelta = sparse(1,1);

sol.dLambda1 = (F2(:))';
sol.dLambda3 = sparse(1,(n+m)*(n+m));
sol.Lambda_Alpha = -1;
sol.Lambda_Beta = sparse(1,1);
sol.Lambda_tDelta = sparse(1,1);
sol.Lambda_Y = sparse(1,n*n);
sol.Data = Lambda1(:)'*gradient.dF2_dD(n,m,B,T,G,Phi);

sol.G3_row_without_Data = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta, sol.Lambda_Y];
end