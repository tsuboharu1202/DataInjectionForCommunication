% % function dF_dD = dF_dD(n,m,T,L)
% % n1 = 2*n+m;

% % fprintf('L: %dx%d\n', size(L,1), size(L,2));
% % L1 = [L';
% %     sparse(n+2*m,T)];
% % L2 = [sparse(n,T);
% %     L';
% %     sparse(2*m,T)];

% % Ez = [speye(n,n),sparse(n,n+m)];
% % Ex = [sparse(n,n),speye(n),sparse(n,m)];
% % Eu = [sparse(m,n),sparse(m,n),speye(m)];


% % Ex1 = [Ex;
% %     sparse(n+m*2,n1)];
% % Ex2 = [sparse(n,n1);
% %     Ex;
% %     sparse(m+m,n1)];
% % Ez2 = [sparse(n,n1);
% %     Ez;
% %     sparse(m+m,n1)];
% % Eu4 = [sparse(2*n+m,n1);
% %     Eu];

% % term1 = kron(L1,Ex1);
% % term2 = kron(L2,Ex2);
% % term3 = kron(L1,Ez2);
% % term4 = kron(Ez2,L1)*implicit.helper.commutation(n1,T);
% % term5 = kron(L1,Eu4);
% % term6 = kron(Eu4,L1)*implicit.helper.commutation(n1,T);

% % dF_dD = term1 + term2 + term3 + term4 + term5 + term6;
% % end

% function dF_dD = dF_dD(n,m,T,L)
% nF  = 2*n + 2*m;

% % --- commutation sanity ---
% Knn = implicit.helper.commutation(n,n);     % vec(A') = Knn*vec(A) for A(nxn)
% Kmn = implicit.helper.commutation(m,n);     % A is (m x n) -> A' is (n x m)

% % --- selectors: Z = Ez*D, X = Ex*D, U = Eu*D ---
% Ez = [speye(n), sparse(n,n+m)];
% Ex = [sparse(n,n), speye(n), sparse(n,m)];
% Eu = [sparse(m,2*n), speye(m)];

% IT = speye(T);
% Sz = kron(IT, Ez);   % vec(Z) = Sz * vec(D)
% Sx = kron(IT, Ex);   % vec(X) = Sx * vec(D)
% Su = kron(IT, Eu);   % vec(U) = Su * vec(D)

% % --- base Jacobians wrt vec(D) ---
% J_ZL  = kron(L', speye(n)) * Sz;     % vec(ZL)
% J_XL  = kron(L', speye(n)) * Sx;     % vec(XL)
% J_UL  = kron(L', speye(m)) * Su;     % vec(UL)

% J_ZLt = Knn * J_ZL;                  % vec((ZL)')
% J_ULt = Kmn * J_UL;                  % vec((UL)')  ★ここが p,q 逆だと壊れる
% J_S   = 0.5*(speye(n*n)+Knn) * J_XL; % vec(S)

% % --- embedding into vec(F1) ---
% Er1 = [speye(n); sparse(nF-n, n)];
% Er2 = [sparse(n,n); speye(n); sparse(2*m, n)];
% Er4 = [sparse(2*n+m, m); speye(m)];

% Ec1 = Er1; Ec2 = Er2; Ec4 = Er4;

% E11 = kron(Ec1, Er1);   % (1,1) block
% E22 = kron(Ec2, Er2);   % (2,2) block
% E21 = kron(Ec1, Er2);   % (2,1) block
% E12 = kron(Ec2, Er1);   % (1,2) block
% E41 = kron(Ec1, Er4);   % (4,1) block
% E14 = kron(Ec4, Er1);   % (1,4) block

% dF_dD = E11*J_S + E22*J_S + E21*J_ZL + E12*J_ZLt + E41*J_UL + E14*J_ULt;
% dF_dD = sparse(dF_dD);
% end
