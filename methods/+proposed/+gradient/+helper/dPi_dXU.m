function [dPi_dX,dPi_dU] = dPi_dXU(n,m,T,Gamma)
tempX = [sparse(m,n);speye(n)];
tempU = [speye(m);sparse(n,m)];
dGamma_dX = kron(speye(T),tempX);
dGamma_dU = kron(speye(T),tempU);


pinvGamma = pinv(Gamma);
dpinvGamma_dGamma = helper.dpinvGamma_dGamma(n,m,T,Gamma);


termX1 = -kron(Gamma',speye(T))*dpinvGamma_dGamma*dGamma_dX;
termX2 = -kron(speye(T),pinvGamma)*dGamma_dX;
dPi_dX = termX1 + termX2;

termU1 = -kron(Gamma',speye(T))*dpinvGamma_dGamma*dGamma_dU;
termU2 = -kron(speye(T),pinvGamma)*dGamma_dU;
dPi_dU = termU1 + termU2;

end