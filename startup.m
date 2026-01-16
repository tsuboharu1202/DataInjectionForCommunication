% startup.m  (ForCommunication 直下)
projroot = fileparts(mfilename('fullpath'));

% ルートを追加（+cfg が見える）
addpath(projroot);

% === 新しいフォルダ構造 ===
addpath(genpath(fullfile(projroot,'core')));
addpath(genpath(fullfile(projroot,'methods')));
addpath(genpath(fullfile(projroot,'attack')));
addpath(genpath(fullfile(projroot,'demos')));
addpath(genpath(fullfile(projroot,'Experiment')));

% === 旧フォルダ構造（アーカイブ、互換性のため一時的に残す） ===
addpath(genpath(fullfile(projroot,'_archive','basic_src')));
addpath(genpath(fullfile(projroot,'_archive','com_src')));
addpath(genpath(fullfile(projroot,'_archive','scripts')));
addpath(genpath(fullfile(projroot,'resources')));

% 衝突回避
clear cfg
clear classes
rehash toolboxcache

% 動作確認（コメントアウト可）
% assert(exist('cfg.Const','class')==8, 'cfg.Const not visible. Check paths.');
