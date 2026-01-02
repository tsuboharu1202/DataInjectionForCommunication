function dF_dD = dF_dD(n,m,T,L)
n1 = 2*n+m;

fprintf('L: %dx%d\n', size(L,1), size(L,2));
L1 = [L';
    sparse(n+2*m,T)];
L2 = [sparse(n,T);
    L';
    sparse(2*m,T)];

Ez = [speye(n,n),sparse(n,n+m)];
Ex = [sparse(n,n),speye(n),sparse(n,m)];
Eu = [sparse(m,n),sparse(m,n),speye(m)];


Ex1 = [Ex;
    sparse(n+m*2,n1)];
Ex2 = [sparse(n,n1);
    Ex;
    sparse(m+m,n1)];
Ez2 = [sparse(n,n1);
    Ez;
    sparse(m+m,n1)];
Eu4 = [sparse(2*n+m,n1);
    Eu];

term1 = kron(L1,Ex1);
term2 = kron(L2,Ex2);
term3 = kron(L1,Ez2);
term4 = kron(Ez2,L1)*implicit.helper.commutation(n1,T);
term5 = kron(L1,Eu4);
term6 = kron(Eu4,L1)*implicit.helper.commutation(n1,T);

dF_dD = term1 + term2 + term3 + term4 + term5 + term6;
end