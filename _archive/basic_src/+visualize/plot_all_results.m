% function plot_all_results(X, dx, U, du, Z, dz)
% % plot_all_results
% %   X, dx, U, du, Z, dz を受け取り、
% %   それぞれ (元) vs (元+差分) の比較グラフを
% %   3つの独立したFigureウィンドウに描画します。

% % --- 1. 状態 (State) の描画 ---
% % 足し合わせた変数を先に作っておくと、plot_data内の凡例(Legend)が
% % "X" と "X_new" のように綺麗に表示されます。
% X_new = X + dx;

% figure('Name', 'State Comparison (X vs X+dx)', 'NumberTitle', 'off');
% visualize.plot_data(X, X_new);

% % --- 2. 入力 (Input) の描画 ---
% U_new = U + du;

% figure('Name', 'Input Comparison (U vs U+du)', 'NumberTitle', 'off');
% visualize.plot_data(U, U_new);

% % --- 3. 出力/評価 (Output) の描画 ---
% Z_new = Z + dz;

% figure('Name', 'Output Comparison (Z vs Z+dz)', 'NumberTitle', 'off');
% visualize.plot_data(Z, Z_new);

% % (オプション) ウィンドウを整列させるならここに配置コードを書く
% end