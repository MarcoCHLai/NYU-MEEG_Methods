close all; clear; clc;

script_dir = '/scratch/cl5564/MEEG_methods/script_ft';
cd(script_dir)

% Initialize Fieldtrip
addpath('/scratch/cl5564/MATLAB/fieldtrip-20240110');
ft_defaults;

%% define directory
data_type = 'EEG';
subj = 'sub-004';

% data directory
data_dir = sprintf('../dataset/%s/%s/%s', data_type, subj, lower(data_type));
cd(data_dir)
data_dir = pwd;

% listing the file names in the folder
flist = dir(data_dir);
flist = {flist.name};

% screening the files based on their extensions
bdf_fname  = fullfile(data_dir, char(flist(endsWith(flist, '.bdf'))));
elec_fname = fullfile(data_dir, char(flist(endsWith(flist, '_electrodes.tsv'))));
eve_fname  = fullfile(data_dir, char(flist(endsWith(flist, '_events.tsv'))));

%% Load data
cfg = [];
cfg.dataset    = bdf_fname;
cfg.continuous = 'yes';
data_raw       = ft_preprocessing(cfg);

%% data segment
% 'sub-004_sess-oddball01.mat': 0-917
% 'sub-004_sess-story01.mat': 921-1600
% 'sub-004_sess-oddball02.mat': 1700-3110
% 'sub-004_sess-story02.mat': 3110-3800
% 'sub-004_sess-oddball03.mat': 3800-4570
cfg = [];
cfg.toilim = [3110 3800];
data = ft_redefinetrial(cfg, data_raw);
save('sub-004_sess-story02.mat', 'data')
