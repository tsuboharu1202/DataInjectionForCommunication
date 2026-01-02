function [dPi_dX,dPi_dU] = dPi_dXU(n,m,T,Gamma)
tempX = [sparse(m,n);speye(n)];
tempU = [speye(m);sparse(n,m)];
dGamma_dX = kron(speye(T),tempX);
dGamma_dU = kron(speye(T),tempU);


pinvGamma = pinv(Gamma);
dpinvGamma_dGamma = implicit_regularization.helper.dpinvGamma_dGamma(n,m,T,Gamma);

% デバッグ用: 各変数のサイズを表示
fprintf('=== dPi_dXU デバッグ情報 ===\n');
fprintf('n=%d, m=%d, T=%d\n', n, m, T);
fprintf('Gamma: %dx%d\n', size(Gamma,1), size(Gamma,2));
fprintf('pinvGamma: %dx%d\n', size(pinvGamma,1), size(pinvGamma,2));
fprintf('dpinvGamma_dGamma: %dx%d\n', size(dpinvGamma_dGamma,1), size(dpinvGamma_dGamma,2));
fprintf('dGamma_dX: %dx%d\n', size(dGamma_dX,1), size(dGamma_dX,2));
fprintf('dGamma_dU: %dx%d\n', size(dGamma_dU,1), size(dGamma_dU,2));
fprintf('kron(Gamma,speye(T)): %dx%d\n', size(kron(Gamma,speye(T)),1), size(kron(Gamma,speye(T)),2));
fprintf('kron(speye(T),pinvGamma): %dx%d\n', size(kron(speye(T),pinvGamma),1), size(kron(speye(T),pinvGamma),2));
fprintf('===========================\n');

termX1 = -kron(Gamma',speye(T))*dpinvGamma_dGamma*dGamma_dX;
termX2 = -kron(speye(T),pinvGamma)*dGamma_dX;
dPi_dX = termX1 + termX2;

termU1 = -kron(Gamma',speye(T))*dpinvGamma_dGamma*dGamma_dU;
termU2 = -kron(speye(T),pinvGamma)*dGamma_dU;
dPi_dU = termU1 + termU2;

end