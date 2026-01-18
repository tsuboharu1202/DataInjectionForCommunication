% function plot_data(varargin)
% % plot_data(X1, X2, ..., XK)
% %  - すべて n×T 行列（同サイズ）を想定
% %  - 各行ごとに K 本の系列を重ね描き
% %  - 凡例のラベルは呼び出し側の「変数名」をそのまま使用
% %
% % 例: plot_data(X_ori, X_ev, X_score)

%     % --- 入力チェック ---
%     assert(nargin>=1, '少なくとも1個の行列を渡してください。');
%     % サイズ確認
%     sz1 = size(varargin{1});
%     for k = 1:nargin
%         assert(ismatrix(varargin{k}) && isnumeric(varargin{k}), ...
%             '第%d引数は数値2次元配列である必要があります。', k);
%         assert(isequal(size(varargin{k}), sz1), ...
%             '全行列は同じサイズ (n×T) である必要があります。');
%     end
%     [n,T] = deal(sz1(1), sz1(2));
%     t = 1:T;

%     % ラベル名（変数名）。式や一時値だと inputname は空→"X1"等でフォールバック
%     labels = cell(1,nargin);
%     for k = 1:nargin
%         nm = inputname(k);
%         if isempty(nm), nm = sprintf('X%d',k); end
%         labels{k} = nm;
%     end

%     % --- 描画 ---
%     tl = tiledlayout(n,1,'TileSpacing','compact','Padding','compact');
%     title(tl,'Time series per state (rows)');

%     styles = {'-','--',':','-.', ...
%               '-','--',':','-.'};    % 足りなければ循環
%     legAx = [];

%     for i = 1:n
%         ax = nexttile;
%         hold(ax,'on');
%         for k = 1:nargin
%             plot(t, varargin{k}(i,:), styles{1+mod(k-1,numel(styles))}, ...
%                  'LineWidth', 1.2, 'DisplayName', labels{k});
%         end
%         grid(ax,'on');
%         ylabel(ax, sprintf('row %d', i));
%         if i==n, xlabel(ax,'k'); end
%         legAx = ax; % 最後の軸を凡例元に
%     end

%     lg = legend(legAx, 'show', 'Orientation','horizontal');
%     try, lg.Layout.Tile = 'north'; end
% end
