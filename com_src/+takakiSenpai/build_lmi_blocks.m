% function [M25, M26, M27] = build_lmi_blocks(Y,L,alpha,beta,tTau,tDelta,H,E1,E2,G,Phi)
%     % sizes
%     n  = size(Y,1);
%     m  = size(L,1);
%     n2 = size(E2,1);

%     EY  = E1*Y + E2*L;
%     HtH = H*H.';            % n x n

%     % shorthand zeros
%     Znn  = zeros(n,n);   Znm  = zeros(n,m);  Zmn  = zeros(m,n);
%     Zmm  = zeros(m,m);   Znn2 = zeros(n,n2); Zn2n = zeros(n2,n);
%     Zmn2 = zeros(m,n2);  Zn2m = zeros(n2,m); Zn2n2 = zeros(n2,n2);

%     % ---------- (25) Left big symmetric block ----------
%     row1 = [ Y - tTau*HtH - beta*eye(n),  Znn,  Znm,    Znm,   Znn,   Znn2,  Znm ];
%     row2 = [ Znn,                          Znn,  Znm,    Znm,  -Y,     Znn2,  Znm ];
%     row3 = [ Zmn,                          Zmn,  Zmm,   -eye(m), -L,    Zmn2,  Zmm ];
%     row4 = [ Zmn,                          Zmn, -eye(m),  eye(m), Zmn,  E2',   Zmm ];
%     row5 = [ Znn,                         -Y',  -L',     Znm,    Y,     EY',   L'   ];
%     row6 = [ Zn2n,                         Zn2n, Zn2m,    E2,     EY,    tTau*eye(n2), zeros(n2,m) ];
%     row7 = [ Zmn,                          Zmn,  Zmm,     Zmm,    L,    zeros(m,n2),  tDelta*eye(m) ];

%     Left = [ row1;
%              row2;
%              row3;
%              row4;
%              row5;
%              row6;
%              row7 ];

%     % data-dependent term: G * Phi * G'
%     M25 = Left - alpha*(G*Phi*G');

%     % ---------- (26) ----------
%     M26 = [ eye(m)    ,     zeros(m,n)  ,    E2'     ,      zeros(m,m);
%             zeros(n,m) ,    Y       ,        EY'        ,     L';
%             E2      ,     EY         ,     tTau*eye(n2)  ,  zeros(n2,m);
%             zeros(m,m)   ,  L         ,      zeros(m,n2)  ,   tDelta*eye(m) ];

%     % ---------- (27) ----------
%     M27 = [ Y            EY'           L';
%             EY           tTau*eye(n2)  zeros(n2,m);
%             L            zeros(m,n2)   tDelta*eye(m) ];
% end
