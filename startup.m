% startup.m  (ForCommunication 直下)
projroot = fileparts(mfilename('fullpath'));

% ルートを追加（+cfg が見える）
addpath(projroot);

% === 外部ツール（ローカル用） ===
tools_dir = '/Users/tsuboharu1202/Documents/MATLAB/_tools';
if exist(tools_dir, 'dir')
    addpath(genpath(fullfile(tools_dir, 'YALMIP')));
    addpath(genpath(fullfile(tools_dir, 'MOSEK', '11.0', 'toolbox')));
    addpath(genpath(fullfile(tools_dir, 'SDPT3')));
    addpath(genpath(fullfile(tools_dir, 'SeDuMi')));
end

% === 新しいフォルダ構造 ===
addpath(genpath(fullfile(projroot,'core')));
addpath(genpath(fullfile(projroot,'methods')));
addpath(genpath(fullfile(projroot,'attack')));
addpath(genpath(fullfile(projroot,'demos')));
addpath(genpath(fullfile(projroot,'Experiment')));

% === 旧フォルダ構造（削除済み - 新しい構造を使用） ===
% addpath(genpath(fullfile(projroot,'_archive','basic_src')));
% addpath(genpath(fullfile(projroot,'_archive','com_src')));
% addpath(genpath(fullfile(projroot,'_archive','scripts')));
addpath(genpath(fullfile(projroot,'resources')));

% 衝突回避
clear cfg
clear classes
rehash toolboxcache

% 動作確認（コメントアウト可）
% assert(exist('cfg.System','class')==8, 'cfg.System not visible. Check paths.');
