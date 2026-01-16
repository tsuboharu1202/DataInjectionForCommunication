function sol = G2_merge(n,m,T,X,Z,U,Pi,gamma,Gamma,L,Lambda,Lambda_P,Lambda3)
% G2 は Lagrangian をLで変微分したもの
sol = struct();
sol.dL = 2*gamma*kron(speye(n),Pi);
dF_dL = implicit_regularization.helper.dF_dL(n,m,T,X,Z,U);
dF_dLambda = -dF_dL';

sol.dLambda = dF_dLambda;
sol.dtDelta = sparse(T*n,1);

sol.dLambda_P = -kron(speye(n),X');
% 後で
dPi_dD = implicit_regularization.helper.dPi_dD(n,m,T,Gamma);
term1 = 2*gamma*kron(L',speye(T))*dPi_dD;
term2 = -implicit_regularization.helper.dVec_Lambda_TdFdL_dD(n,m,T,Lambda);

C2nmT = implicit.helper.commutation(2*n+m,T);
Ex = [sparse(n,n), speye(n),sparse(n,m)];
KTx = implicit.helper.commutation(T,n);
% term3 = - (kron(Lambda_P'*Ex,speye(T)))*C2nmT;
term3_Data = -kron(Lambda_P', speye(T)) * KTx * kron(speye(T), Ex);


Knn = implicit.helper.commutation(n,n);
sol.dLambda3 = kron(speye(n), X') * (speye(n*n) - Knn);

sol.Data = term1 + term2 + term3_Data;
% --- new: Data term due to X dependence of X'*(Lambda3-Lambda3') ---
M   = Lambda3 - Lambda3';
KTx = implicit.helper.commutation(T,n);        % vec(dX') = KTx * vec(dX)
Ex  = [sparse(n,n), speye(n), sparse(n,m)];    % picks X from D=[Z;X;U]
extra_Data = kron(M', speye(T)) * KTx * kron(speye(T), Ex);

sol.Data = sol.Data + extra_Data;

sol.G2_row_without_Data = [sol.dtDelta, sol.dL, sol.dLambda,sol.dLambda_P,sol.dLambda3];

end