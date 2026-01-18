% function sol = G5_merge(n,m,T,B)
% % G5 は Lagrangian をtDeltaで変微分したもの
% sol = struct();
% sol.dL = sparse(1,n*m);
% sol.dY = sparse(1,n*n);
% sol.dAlpha = sparse(1,1);
% sol.dBeta = sparse(1,1);
% sol.dtDelta = sparse(1,1);


% E = [B;
%     zeros(2*n+m,m)];
% M = E*E';
% sol.dLambda1 = M(:)';
% sol.dLambda3 = sparse(1,(n+m)*(n+m));
% sol.Lambda_Alpha = sparse(1,1);
% sol.Lambda_Beta = sparse(1,1);
% sol.Lambda_tDelta = -1;
% sol.Lambda_Y = sparse(1,n*n);
% sol.Data = sparse(1,T*(2*n+m));

% sol.G5_row_without_Data = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta, sol.Lambda_Y];
% end