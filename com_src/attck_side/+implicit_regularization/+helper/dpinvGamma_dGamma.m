function dpinvGamma_dGamma = dpinvGamma_dGamma(n,m,T,Gamma)

Cnm_T = implicit.helper.commutation(n+m,T);
GGT = Gamma*Gamma';
GGT_inv = GGT^(-1);

term1 = kron(GGT_inv',speye(T))*Cnm_T;
term2 = -kron(speye(n+m),Gamma')*kron(GGT_inv',GGT_inv)*(kron(Gamma,speye(m+n)) + kron(speye(m+n),Gamma))*Cnm_T;
dpinvGamma_dGamma = term1 + term2;
end