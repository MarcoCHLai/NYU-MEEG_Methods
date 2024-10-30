close all; clear; clc;

script_dir = '/scratch/cl5564/MEEG_methods/script_ft';
cd(script_dir)

% Initialize Fieldtrip
addpath('/scratch/cl5564/MATLAB/fieldtrip-20240110');
ft_defaults;

%% define directory
data_type = 'EEG';
subj = 'sub-004';
sess = 'sess-oddball';

% data directory
data_dir = sprintf('../dataset/%s/%s/%s', data_type, subj, lower(data_type));
cd(data_dir)
data_dir = pwd;

% listing the file names in the folder
flist = dir(data_dir);
flist = {flist.name};

% load epoch data
fname = fullfile(data_dir, sprintf('%s_%s_data_pruned.mat', subj, sess));
load(fname)

% load event info
eve_fname  = fullfile(data_dir, char(flist(endsWith(flist, '_events.tsv'))));
eve_tab = readtable(eve_fname, 'FileType', 'text');

% read in electrodes with known coordinates
elec_fname = fullfile(data_dir, char(flist(endsWith(flist, '_electrodes.tsv'))));
elec_tab   = readtable(elec_fname, 'FileType', 'text');
elec_tab   = rmmissing(elec_tab); % Remove NaN (for EXG channels)

%% prepare layout and channel neighbor
% prepare electrode info
elec = [];
elec.elecpos = [elec_tab.x, elec_tab.y, elec_tab.z];
elec.chanpos = elec.elecpos;
elec.unit    = 'mm';
elec.label   = data_pruned.label;

% prepare layout
cfg = [];
cfg.layout   = 'biosemi256';
cfg.feedback = 'no';
layout       = ft_prepare_layout(cfg);

% prepare channel neighbor
cfg = [];
cfg.method   = 'triangulation';
cfg.feedback = 'no';
cfg.elec     = elec;
neighbours   = ft_prepare_neighbours(cfg);

%% update trialinfo
% add stimulus category and name

for i = 1:length(data_pruned.trialinfo)
    
    block_num = data_pruned.trialinfo{i}.blocknumber;
    trial_num = data_pruned.trialinfo{i}.trialnumber;
    event_idx = find((eve_tab.blockNumber == block_num) & (eve_tab.trialNumber == trial_num));
    stim_cate = eve_tab(event_idx,:).stimulusCategory;
    stim_name = eve_tab(event_idx,:).stimulusName;
    data_pruned.trialinfo{i}.stimCategory = stim_cate;
    data_pruned.trialinfo{i}.stimName = stim_name;
    
end

%% processing for ERP
cfg = [];
cfg.lpfilter = 'yes';
cfg.lpfreq   = 15;
cfg.detrend  = 'yes';
cfg.demean   = 'yes';           % baseline correcttion
cfg.baselinewindow = [-0.5 0];  % using the mean activity in this window
data_bl = ft_preprocessing(cfg, data_pruned);

%%
stim_type = 'classicalAud'; % classicialAud, semanticAud, semanticVis

% odd condtion
cfg = [];
cfg.trials = find(cellfun(@(x) ...
             strcmp(stim_type, string(x.blockType)) && ...
             strcmp('odd', string(x.type)), ...
             data_bl.trialinfo));
cfg.keeptrials = 'yes';
evk_odd = ft_timelockanalysis(cfg, data_bl);

% even condtion
cfg = [];
cfg.trials = find(cellfun(@(x) ...
             strcmp(stim_type, string(x.blockType)) && ...
             strcmp('even', string(x.type)), ...
             data_bl.trialinfo));
cfg.keeptrials = 'yes';
evk_even = ft_timelockanalysis(cfg, data_bl);

