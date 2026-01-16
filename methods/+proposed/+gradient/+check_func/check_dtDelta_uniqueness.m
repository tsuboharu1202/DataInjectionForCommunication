function check_dtDelta_uniqueness(H, n, T, n1)
% 変数の並び: [tDelta(1), vecL(T*n), vecLambda(n1^2), vecLambdaP(n^2), vecLambda3(n^2)]
s = svd(full(H),'econ');
tol = 1e-12*s(1);

[U,S,V] = svd(full(H),'econ');
r = sum(diag(S) > tol);
if r == size(H,2)
    fprintf("Full column rank: unique solution.\n");
    return;
end

Vnull = V(:, r+1:end); % basis
fprintf("null-dim = %d\n", size(Vnull,2));

% dtDelta 成分
dt = Vnull(1,:);
fprintf("max |null(dtDelta)| = %.3e\n", max(abs(dt)));

% ブロックごとの寄与を見る（どこが自由度源か）
idx = 1;
idx_t = idx; idx = idx+1;
idx_L = idx:(idx+T*n-1); idx = idx+T*n;
idx_Lam = idx:(idx+n1*n1-1); idx = idx+n1*n1;
idx_LamP = idx:(idx+n*n-1); idx = idx+n*n;
idx_Lam3 = idx:(idx+n*n-1);

fprintf("||null(L)||=%.3e\n", norm(Vnull(idx_L,:), 'fro'));
fprintf("||null(Lambda)||=%.3e\n", norm(Vnull(idx_Lam,:), 'fro'));
fprintf("||null(LambdaP)||=%.3e\n", norm(Vnull(idx_LamP,:), 'fro'));
fprintf("||null(Lambda3)||=%.3e\n", norm(Vnull(idx_Lam3,:), 'fro'));
end
