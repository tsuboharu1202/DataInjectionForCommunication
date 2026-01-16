function dVec_Lambda_TdFdL_dD = dVec_Lambda_TdFdL_dD(n,m,T,Lambda_ori)

Lambda_sym = 0.5*(Lambda_ori + Lambda_ori');
n1 = 2*n+m;
Lambda11 = Lambda_sym(1:n,1:n);
Lambda21 = Lambda_sym(n+1:n+n,1:n);
Lambda22 = Lambda_sym(n+1:n+n,n+1:n+n);
Lambda41 = Lambda_sym(2*n+m+1:2*n+2*m,1:n);


Ez = [speye(n,n),sparse(n,n+m)];
Ex = [sparse(n,n),speye(n),sparse(n,m)];
Eu = [sparse(m,n),sparse(m,n),speye(m)];

LambdaX = (Lambda11 + Lambda22)'*Ex;
LambdaZ = (Lambda21')*Ez;
LambdaU = (Lambda41')*Eu;

I_T = speye(T);

Cn1T = implicit.helper.commutation(n1,T);

term1 = kron(LambdaX,I_T)*Cn1T;
term2 = 2*kron(LambdaZ,I_T)*Cn1T;
term3 = 2*kron(LambdaU,I_T)*Cn1T;

dVec_Lambda_TdFdL_dD = term1 + term2 + term3;

end