function dF1_dL = dF1_dL(n,m, B)
E1 = [speye(n);
    zeros(2*n+m,n)];

E2 = [zeros(2*n,n);
    speye(n);
    zeros(m,n)];

E3 = [zeros(3*n,m);
    speye(m)];


term1 = kron(E1,E2*B);
term2 = kron(E2*B,E1)*implicit.helper.commutation(m,n);
term3 = kron(E3,E2);
term4 = kron(E2,E3)*implicit.helper.commutation(n,m);

dF1_dL = term1 + term2 + term3 + term4;