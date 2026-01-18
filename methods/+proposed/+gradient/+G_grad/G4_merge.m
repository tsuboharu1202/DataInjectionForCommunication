function sol = G4_merge(n,m,T)
% G4 は Lambda - Lambda' = 0 をLで変微分したもの
n1 = 2*n+2*m;
n1_2 = n1*n1;
sol = struct();
sol.dL = sparse(n1_2,n*T);
sol.dLambda = speye(n1_2) - core.helper.commutation(n1,n1);
sol.dtDelta = sparse(n1_2,1);
sol.dLambda_P = sparse(n1_2,n*n);
sol.Data = sparse(n1_2,T*(2*n+m));
sol.dLambda3 = sparse(n1*n1, n*n);
sol.G4_row_without_Data = [sol.dtDelta, sol.dL, sol.dLambda, sol.dLambda_P, sol.dLambda3];


end