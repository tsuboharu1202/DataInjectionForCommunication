function dF1_dL = dF1_dL(n,m, B)
E1 = [speye(n);
    zeros(2*n+m,n)];

E3= [zeros(2*n,n);
    speye(n);
    zeros(m,n)];

E4 = [zeros(3*n,m);
    speye(m)];


term1 = kron(E3,E1*B);
term2 = kron(E1*B,E3)*core.helper.commutation(m,n);
term3 = kron(E3,E4);
term4 = kron(E4,E3)*core.helper.commutation(m,n);

dF1_dL = term1 + term2 + term3 + term4;