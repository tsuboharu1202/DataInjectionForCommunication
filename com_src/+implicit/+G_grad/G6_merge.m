function sol = G6_merge(n,m,T,Lambda1,alpha,F1,F2,B,G,Phi)
% G6 は 相補性条件における Lambda1*(F1-alpha*F2) = 0
% 注意: F2 = G*Phi*G' (alphaを含まない) なので、F1 - alpha*F2 = F1 - alpha*(G*Phi*G')
n1_2 = (3*n+m)*(3*n+m);
n2_2 = (n+m)*(n+m);
eye = speye(3*n+m);
sub_matrix = kron(Lambda1,eye);

sol = struct();
sol.dL = sub_matrix*implicit.helper.dF1_dL(n,m,B);
sol.dY = sub_matrix*implicit.helper.dF1_dY(n,m);

% d/dalpha [Lambda1 * (F1 - alpha*F2)] = Lambda1 * d/dalpha [F1 - alpha*F2]
% d/dalpha [F1 - alpha*F2] = -F2
% したがって、d/dalpha [Lambda1 * (F1 - alpha*F2)] = -Lambda1 * F2
% vec化すると: d/dalpha [vec(Lambda1 * (F1 - alpha*F2))] = -vec(Lambda1 * F2)
F2Lambda1 = F2*Lambda1;
sol.dAlpha = -F2Lambda1(:);

E11 = [speye(n);
    zeros(2*n+m,n)];
M11 = E11*E11';
M11Lambda1 = -M11*Lambda1;

sol.dBeta = M11Lambda1(:);

E11_withB = [B;
    zeros(2*n+m,m)];
M11_withB = E11_withB*E11_withB';
M11_withBLambda1 = -M11_withB*Lambda1;

sol.dtDelta = M11_withBLambda1(:);


% 修正: d/dLambda1 [Lambda1 * (F1 - alpha*F2)] = (F1 - alpha*F2) をvec化
% G14と同様に、kron(speye, F1_minus_alphaF2)*commutation を使う
F1_minus_alphaF2 = F1 - alpha*F2;  % F2 = G*Phi*G' (alphaを含まない) なので、F1 - alpha*F2 = F1 - alpha*(G*Phi*G')
% vec(Lambda1 * M) = (M^T \otimes I) vec(Lambda1) = K (I \otimes M) vec(Lambda1)
% G14の形式に合わせる: kron(speye(n), Y) * commutation
% d/dvec(Lambda1) [vec(Lambda1 * (F1 - alpha*F2))] = kron(speye, F1_minus_alphaF2) * commutation
sol.dLambda1 = kron(speye(3*n+m), F1_minus_alphaF2)*implicit.helper.commutation(3*n+m,3*n+m);
sol.dLambda3 = sparse(n1_2,n2_2);
sol.Lambda_Alpha = sparse(n1_2,1);
sol.Lambda_Beta = sparse(n1_2,1);
sol.Lambda_tDelta = sparse(n1_2,1);
sol.Lambda_Y = sparse(n1_2,n*n);
sol.Data = -alpha*sub_matrix*implicit.helper.dF2_dD(n,m,B,T,G,Phi);

sol.G6_row_without_Data = [sol.dL, sol.dY, sol.dAlpha, sol.dBeta, sol.dtDelta, sol.dLambda1, sol.dLambda3, sol.Lambda_Alpha, sol.Lambda_Beta, sol.Lambda_tDelta, sol.Lambda_Y];
end