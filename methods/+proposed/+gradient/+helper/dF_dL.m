% function dF_dL = dF_dL(n,m,T,X,Z,U)
% E1 = [speye(n);
%     sparse(n+2*m,n)];
% E2 = [sparse(n,n);
%     speye(n);
%     sparse(2*m,n)];
% X1 = [X;
%     sparse(2*m+n,T)];
% X2 = [sparse(n,T);
%     X;
%     sparse(2*m,T)];
% Z2 = [sparse(n,T);
%     Z;
%     sparse(2*m,T)];
% U4 = [sparse(2*n+m,T);
%     U];
% dF_dL = kron(E1,X1)+ kron(E2,X2) + kron(E2,Z2) +...
%     kron(E1,Z2) +kron(Z2,E1)*core.helper.commutation(T,n) +...
%     kron(E1,U4) +kron(U4,E1)*core.helper.commutation(T,n);

% end


function dF_dL = dF_dL(n,m,T,X,Z,U)
nF = 2*n + 2*m;

Knn = core.helper.commutation(n,n);
Kmn = core.helper.commutation(m,n);  % ★修正: (m,n)

J_ZL  = kron(speye(n), Z);        % vec(ZL)  = (I ⊗ Z) vec(L)
J_XL  = kron(speye(n), X);        % vec(XL)  = (I ⊗ X) vec(L)
J_UL  = kron(speye(n), U);        % vec(UL)  = (I ⊗ U) vec(L), UL is m×n

J_ZLt = Knn * J_ZL;               % vec((ZL)') = K_{n,n} vec(ZL)
J_ULt = Kmn * J_UL;               % ★ vec((UL)') = K_{m,n} vec(UL)
J_S   = 0.5*(speye(n*n)+Knn)*J_XL;

% row selectors
Er1 = [speye(n);           sparse(nF-n, n)];
Er2 = [sparse(n,n);        speye(n); sparse(2*m, n)];
Er4 = [sparse(2*n+m, m);   speye(m)];

Ec1 = Er1; Ec2 = Er2; Ec4 = Er4;

E11 = kron(Ec1, Er1);   % (1,1)
E22 = kron(Ec2, Er2);   % (2,2)
E21 = kron(Ec1, Er2);   % (2,1)
E12 = kron(Ec2, Er1);   % (1,2)
E41 = kron(Ec1, Er4);   % (4,1)
E14 = kron(Ec4, Er1);   % (1,4)

dF_dL = E11*J_S + E22*J_S + E21*J_ZL + E12*J_ZLt + E41*J_UL + E14*J_ULt;
dF_dL = sparse(dF_dL);
end
