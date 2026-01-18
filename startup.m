% startup.m  (ForCommunication 直下)
projroot = fileparts(mfilename('fullpath'));

% ルートを追加（+cfg が見える）
addpath(projroot);


% --- tools_dir を local / online で自動切替 ---
% --- tools_dir を local / online で自動切替（堅牢版） ---
tools_dir_local = '/Users/tsuboharu1202/Documents/MATLAB/_tools';

% Onlineのドライブ根を推定（MATLAB_DRIVEが空でも動く）
matlab_drive = getenv('MATLAB_DRIVE');   % 取れれば使う
if isempty(matlab_drive)
    matlab_drive = userpath;            % ここが /MATLAB Drive になってる
end

% userpath が "pathsep" 区切りで複数返す場合に備える
matlab_drive = strsplit(matlab_drive, pathsep);
matlab_drive = matlab_drive{1};

tools_dir_online = fullfile(matlab_drive, '_tools');

tools_dir = '';
if exist(tools_dir_local, 'dir')
    tools_dir = tools_dir_local;
elseif exist(tools_dir_online, 'dir')
    tools_dir = tools_dir_online;
end


if exist(tools_dir, 'dir')
    addpath(genpath(fullfile(tools_dir, 'YALMIP')));
    addpath(genpath(fullfile(tools_dir, 'MOSEK', '11.0', 'toolbox')));
    addpath(genpath(fullfile(tools_dir, 'SDPT3')));
    addpath(genpath(fullfile(tools_dir, 'SeDuMi')));
end

% === 新しいフォルダ構造 ===
addpath(genpath(fullfile(projroot,'+core')));
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
