% function sol = G1_merge(n,m, B,T)
% % G1 は Lagrangian をdeltaで変微分したもの
% sol = struct();
% sol.dL = sparse(1,n*T);
% dF_dDelta = [sparse(n,n), sparse(n,n), sparse(n,m), sparse(n,m);
%     sparse(n,n), sparse(n,n), B, sparse(n,m);
%     sparse(m,n), B', sparse(m,m), sparse(m,m);
%     sparse(m,n), sparse(m,n), sparse(m,m), sparse(m,m)];

% sol.dLambda = -dF_dDelta(:)';
% sol.dtDelta = 0;
% sol.dLambda_P = sparse(1,n*n);
% sol.Data = sparse(1,(2*n+m)*T);
% sol.dLambda3 = sparse(1, n*n);
% sol.G1_row_without_Data = [sol.dtDelta, sol.dL, sol.dLambda, sol.dLambda_P, sol.dLambda3];
% % sol.Data はそのまま

% end