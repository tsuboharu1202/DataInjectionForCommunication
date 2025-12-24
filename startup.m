% startup.m  (ForCommunication 直下)
projroot = fileparts(mfilename('fullpath'));

% 親だけ通す（+cfg 自体は通さない）
addpath(projroot);
addpath(genpath(fullfile(projroot,'basic_src')));
addpath(genpath(fullfile(projroot,'resources')));
addpath(genpath(fullfile(projroot,'scripts')));

% パッケージ（+で始まるフォルダ）を先に追加して優先順位を上げる
addpath(genpath(fullfile(projroot,'com_src','+attack')));
addpath(genpath(fullfile(projroot,'com_src','+original_thesis')));
% その後、com_src全体を追加（ForTakakiSenpaiThesisなども含む）
addpath(genpath(fullfile(projroot,'com_src')));


% 衝突回避
clear cfg
clear classes
clear functions  % 関数キャッシュをクリア（ForLQRとの混同を防ぐ）
rehash toolboxcache

% 動作確認（コメントアウト可）
% assert(exist('cfg.Const','class')==8, 'cfg.Const not visible. Check paths.');
