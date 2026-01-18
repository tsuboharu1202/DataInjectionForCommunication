function dPi_dD = dPi_dD(n,m,T,Gamma)
Ex = [sparse(n,n),speye(n),sparse(n,m)];
Eu = [sparse(m,n),sparse(m,n),speye(m)];


[dPi_dX,dPi_dU] = proposed.gradient.helper.dPi_dXU(n,m,T,Gamma);
dPi_dD = dPi_dX*kron(speye(T),Ex) + dPi_dU*kron(speye(T),Eu);
end