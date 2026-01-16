function [F1, F2, F3] = build_lmi_blocks(Y,L,alpha,beta,tDelta,G,Phi,sd)
% sizes
n  = size(Y,1);
m  = size(L,1);

B = sd.B;

% shorthand zeros
Znn  = zeros(n,n);   Znm  = zeros(n,m);
Zmn  = zeros(m,n);

% ---------- (12) Left big symmetric block ----------
% 論文の数式(12)に従う:
% - row1の(1,1)ブロック: Y - δ²BB^T - βI, ここでtDelta = δ²なので tDelta*B*B'
% - F2の符号: 論文では F1 - α * G * Phi * G' >= 0
%   したがって、F2 = alpha*(G*Phi*G')として、F1 - F2 >= 0が正しい制約になる
row1 = [ Y - tDelta*(B*B') - beta*eye(n),  Znn,  B*L,    Znm ];
row2 = [ Znn,                          Znn,  Y,    Znm ];
row3 = [ L'*B',                          Y',  Y,    L' ];
row4 = [ Zmn,                          Zmn,  L, eye(m) ];

F1 = [ row1;
    row2;
    row3;
    row4 ];

% data-dependent term: G * Phi * G'
F2 = G*Phi*G';

F3 = [Y , L';
    L, eye(m) ];
end
