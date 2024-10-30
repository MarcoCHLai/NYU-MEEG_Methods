close all; clear; clc;

script_dir = '/scratch/cl5564/MEEG_methods/script_ft';
cd(script_dir)

% Initialize Fieldtrip
addpath('/scratch/cl5564/MATLAB/fieldtrip-20240110');
ft_defaults;

%% define directory
data_type = 'EEG';
subj = 'sub-004';
sess = 'sess-oddball03';

% data directory
data_dir = sprintf('../dataset/%s/%s/%s', data_type, subj, lower(data_type));
cd(data_dir)
data_dir = pwd;

% listing the file names in the folder
flist = dir(data_dir);
flist = {flist.name};

% screening the files based on their extensions
% bdf_fname  = fullfile(data_dir, char(flist(endsWith(flist, '.bdf'))));
elec_fname = fullfile(data_dir, char(flist(endsWith(flist, '_electrodes.tsv'))));
eve_fname  = fullfile(data_dir, char(flist(endsWith(flist, '_events.tsv'))));
data_fname = fullfile(data_dir, sprintf('%s_%s.mat', subj, sess));
load(data_fname)

%% read events
eve_tab = readtable(eve_fname, 'FileType', 'text');

%% read in electrodes with known coordinates
elec_tab = readtable(elec_fname, 'FileType', 'text');
elec_tab = rmmissing(elec_tab); % Remove NaN (for EXG channels)

%% Keep the electrode positions
n_elec   = size(elec_tab,1);
n_sample = size(data.trial{1},2);
SR       = data.fsample;

elec = [];
elec.elecpos = [elec_tab.x, elec_tab.y, elec_tab.z];
elec.chanpos = elec.elecpos;
elec.unit    = 'mm';
elec.label   = data.label(1:size(elec.elecpos,1));
data.elec    = elec;

%% prepare the standard layout (we can do this from elec too but need to adjust scale and rotation)
cfg = [];
cfg.layout   = 'biosemi256';
cfg.feedback = 'yes';
layout       = ft_prepare_layout(cfg);

%% re-reference the data
cfg = [];
cfg.channel    = data.label(1:n_elec); % this is the default
cfg.reref      = 'yes';
cfg.refmethod  = 'avg';
cfg.refchannel = data.label(1:n_elec);
data_reref     = ft_preprocessing(cfg, data);

%% plotting data via ft_databrowser
cfg = [];
cfg.layout           = layout;
cfg.viewmode         = 'vertical';
cfg.preproc.hpfilter = 'yes';
cfg.preproc.hpfreq   = 0.5;
cfg.ylim             = [-20 20];
cfg.blocksize        = 10;
cfg.position         = [100 250 1000 600];
ft_databrowser(cfg, data_reref);

%% identify bad channels
% create a structure array for storing the bad channel names
bad_channels = {'A26','B7','B8','B9','B19',...
                'C9','C30','D15','E5',...
                'F2','F22','G25','G29'};

% bad channels from the original script in the class
%'B7, 9, 25?, C18?, 'G30'? !F2 E5
%{
bad_channels = {'A26', 'B2', 'B3', 'B4', 'B7' 'B8',...
    'B19', 'C9', 'C30','D15','E5', 'F2' ,'F22',...
    'G25', 'G29'};
%}

% go back one step and interpolate
% then re-reference again!
% then epoch the data

%% prepare a neighborhood structure for interpolation
%(we need to know what channels to use)
cfg = [];
cfg.method   = 'triangulation';
cfg.feedback = 'yes';
cfg.elec     = elec;
neighbours   = ft_prepare_neighbours(cfg);

%% Interpolate bad channels
cfg = [];
cfg.method     = 'weighted';
cfg.neighbours = neighbours;
cfg.badchannel = bad_channels; % List your bad channels
data_interp    = ft_channelrepair(cfg, data);

%% redo the average reference
cfg = [];
cfg.channel    = data.label(1:n_elec); 
cfg.reref      = 'yes';
cfg.refmethod  = 'avg';
cfg.refchannel = data.label(1:n_elec);
data_reref     = ft_preprocessing(cfg, data_interp);

% compute the average here
avg_vec = mean(data_interp.trial{1}(1:n_elec,:));
% clean up
clear data_interp;

%% concatenate with the EOG;
cfg = [];
cfg.channel = data.label(n_elec+1:n_elec+4);

data_eye = ft_selectdata(cfg, data);

% because eog channels are recorded with amplifier refrence, we need to
% re-reference them as well
data_eye.trial{1} = data_eye.trial{1} - avg_vec;
% EOG channels are recognized by databrowser via ft_channelselection
data_eye.label = {'EOG-1'; 'EOG-2'; 'EOG-3';'EOG-4'};
data_combined = ft_appenddata([],data_reref, data_eye);

%% now manually mark artifacts
cfg = [];
cfg.layout           = layout;
cfg.viewmode         = 'vertical';
cfg.preproc.hpfilter = 'yes';
cfg.preproc.hpfreq   = 0.5;
cfg.ylim             = [-20 20];
cfg.blocksize        = 10;
cfg.eogscale         = 0.8; % Four channels are eye channels
cfg.position         = [100 250 1000 600];
cfg.artifactalpha    = 0.8;
cfg_art = ft_databrowser(cfg, data_combined);

%{
% save artifacts info
save(sprintf('art_%s.mat', sess), 'cfg_art')
%}

%% compute one big ICA ...
sample_begin = data_combined.sampleinfo(1);
artifact_def = cfg_art.artfctdef.visual.artifact;
art_vec = zeros(1, n_sample);
for aa = 1 : size(artifact_def,1)
    idx1 = artifact_def(aa,1)-sample_begin+1;
    idx2 = artifact_def(aa,2)-sample_begin+1;
    art_vec(idx1:idx2) = 1;
end

% get the rank of the data (is 240 here because we 
% excluded -15 channels and did avg reference -1)

%% ... without eye channels
rnk = rank(data_reref.trial{1});

% highpass filter for ica
cfg = [];
cfg.hpfilter = 'yes';
cfg.hpfreq   = 1;
data_nan     = ft_preprocessing(cfg, data_reref);

% mask artifacts for ICA
data_nan.trial{1}(:,art_vec==1) = nan;

% compute the ICA on masked+filtered data
cfg = [];
cfg.method       = 'runica';
cfg.numcomponent = rnk;
data_comp        = ft_componentanalysis(cfg, data_nan);

%%
cfg = [];
cfg.layout   = layout;
cfg.viewmode = 'component';
ft_databrowser(cfg, data_comp);
% remove 1 5? 7

%% ... or with eye channels
rnk = rank(data_combined.trial{1});

% highpass filter for ica
cfg = [];
cfg.hpfilter = 'yes';
cfg.hpfreq   = 1;
data_nan     = ft_preprocessing(cfg, data_combined);

% mask artifacts for ICA
data_nan.trial{1}(:,art_vec==1) = nan;

% compute the ICA on masked+filtered data
cfg = [];
cfg.method       = 'runica';
cfg.numcomponent = rnk;
cfg.randomseed   = 1119;
data_comp        = ft_componentanalysis(cfg, data_nan);

%{
% save ICA component
save(sprintf('ICAcomp_%s.mat', sess), 'data_comp')
%}

%% add eye channels to layout
% add eye positions to layout
layeye        = layout;
layeye.label  = [layout.label; {'EOG-1'; 'EOG-2'; 'EOG-3'; 'EOG-4'}];
% layeye.pos    = [lay.pos; [-0.4 0.28]; [-0.2  0.4]; [0.2 0.4]; [0.4 0.28]]; % Example positions
layeye.pos    = [layout.pos; [-0.49 0.49]; [-0.48 0.48];  [0.48 0.48]; [0.49 0.49]]; % Example positions
layeye.width  = [layout.width;ones(4, 1).*layout.width(1)];
layeye.height = [layout.width;ones(4, 1).*layout.width(1)];

cfg = [];
cfg.layout = layeye;
ft_layoutplot(cfg)

%% plotting data via ft_databrowser
cfg = [];
cfg.layout    = layeye;
cfg.viewmode  = 'component';
cfg.blocksize = 10;
cfg.position  = [100 250 1000 600];
ft_databrowser(cfg, data_comp);

%% use unmixing matrix on unfiltered data
unmixing_eye = data_comp.unmixing;

cfg = [];
cfg.unmixing  = unmixing_eye;
cfg.topolabel = data_combined.label;
data_comp     = ft_componentanalysis(cfg, data_combined);

%% reject components:
% oddball-sess1: 1 3 4 5 6 13
% oddball-sess2: 1 3 8 12
% oddball-sess2: 1 2 9 12
cfg           = [];
cfg.component = [1 2 9 12];
data_clean    = ft_rejectcomponent(cfg, data_comp, data_combined);

%% plotting data via ft_databrowser
cfg = [];
cfg.layout           = layeye;
cfg.viewmode         = 'vertical';
cfg.preproc.hpfilter = 'yes';
cfg.preproc.hpfreq   = 0.5;
cfg.ylim             = [-20 20];
cfg.blocksize        = 10;
cfg.position         = [100 250 1000 600];
ft_databrowser(cfg, data_clean)

%%
fname = fullfile(data_dir, sprintf('%s_%s_data_clean.mat', subj, sess));
save(fname, 'data_clean', '-v7.3');
clear data_comp data_combined data_reref data

%% organizing event information
% block code
% classicalAud: 1; semanticVis:  2; semanticAud:  4
% condition code
% even: 64; odd: 128

% 'sub-004_sess-oddball01.mat': 0-917
% 'sub-004_sess-story01.mat': 921-1600
% 'sub-004_sess-oddball02.mat': 1700-3110
% 'sub-004_sess-story02.mat': 3110-3800
% 'sub-004_sess-oddball03.mat': 3800-4570

if strcmp(sess, 'sess-oddball01')
    sess_idx = (eve_tab.timeOnset >= 0) & (eve_tab.timeOnset <= 917);
elseif strcmp(sess, 'sess-oddball02')
    sess_idx = (eve_tab.timeOnset >= 1700) & (eve_tab.timeOnset <= 3110);
elseif strcmp(sess, 'sess-oddball03')
    sess_idx = (eve_tab.timeOnset >= 3800) & (eve_tab.timeOnset <= 4570);
end

% event info for a specfic session
sess_eve = eve_tab(sess_idx,:);
n_eve = size(sess_eve,1);
n_trl = sum(sess_eve.trialNumber~=0);
samp  = nan(n_trl, 1);
info  = cell(n_trl,1);
inf.group = '';
inf.type  = '';
trl_idx   = 0;

for ee = 1: n_eve
    switch sess_eve.trigger_codes_fixed(ee)
        case 1
            inf.group = 'tone';
        case 2
            inf.group = 'vis';
        case 4
            inf.group = 'aud';
        case 8
            inf.group = 'story_start';
        case 16
            inf.group = 'story_seg';
        case 64
            trl_idx = trl_idx+1;
            samp(trl_idx) = sess_eve.sampleOnset(ee);
            info{trl_idx} = inf;
            info{trl_idx}.blockType = sess_eve.blockType(ee);
            info{trl_idx}.trialnumber = sess_eve.trialNumber(ee);
            info{trl_idx}.sampleinfo = sess_eve.sampleOnset(ee);
            info{trl_idx}.type  = 'even';
            info{trl_idx}.blocknumber = sess_eve.blockNumber(ee);
            info{trl_idx}.trigger = sess_eve.trigger_codes_fixed(ee);
        case 128
            trl_idx = trl_idx+1;
            samp(trl_idx) = sess_eve.sampleOnset(ee);
            info{trl_idx} = inf;
            info{trl_idx}.blockType = sess_eve.blockType(ee);
            info{trl_idx}.trialnumber = sess_eve.trialNumber(ee);
            info{trl_idx}.sampleinfo = sess_eve.sampleOnset(ee);
            info{trl_idx}.type  = 'odd';
            info{trl_idx}.blocknumber = sess_eve.blockNumber(ee);
            info{trl_idx}.trigger = sess_eve.trigger_codes_fixed(ee);     
    end
end

%% epoch segment
% cut into 3 second trials (lots of ovelap)
% cfg.trl=[begsample endsample offset]
% https://www.fieldtriptoolbox.org/faq/what_is_the_relation_between_events_such_as_triggers_and_trials/
cfg = [];
cfg.trl  = [samp-1*SR, samp+2*SR, repmat(-1*SR, [n_trl, 1])];
data_cut = ft_redefinetrial(cfg, data_clean);
clear data_clean;

%% remove eye channels
cfg = [];
cfg.channel = data_cut.label(1:n_elec);
data_cut = ft_selectdata(cfg, data_cut);

%% inspect one more time
cfg = [];
cfg.preproc.detrend = 'yes'; % now we detrend instead of highpass filtering! 
cfg.viewmode  = 'vertical';
cfg.ylim      = [-10 10];
ft_databrowser(cfg, data_cut);

%% overwrite trialinfo
data_cut.trialinfo = info;

%% manually identify bad epochs from the data browser
cfg = [];
cfg.preproc.detrend = 'yes'; % now we detrend instead of highpass filtering! 
cfg.viewmode        = 'butterfly';
cfg.ylim            = [-1000 1000];
cfg.position        = [100 250 1000 600];
cfg.artifactalpha   = 0.8;
cfg = ft_databrowser(cfg, data_cut);

%% reject the identified bad epochs via ft_rejectartifact
data_pruned = ft_rejectartifact(cfg, data_cut);

%{
% save data
save(sprintf('%s_%s_data_pruned.mat', subj, sess), 'data_pruned')
%}

%% concatenate data
load('sub-004_sess-oddball01_data_pruned.mat')
data_sess1 = data_pruned;
load('sub-004_sess-oddball02_data_pruned.mat')
data_sess2 = data_pruned;
load('sub-004_sess-oddball03_data_pruned.mat')
data_sess3 = data_pruned;

% concatenate 3 oddball sessions together
cfg = [];
data_pruned = ft_appenddata(cfg, data_sess1, data_sess2, data_sess3);

%{
% save data file
save(sprintf('%s_sess-oddball_data_pruned.mat', subj), 'data_pruned', '-v7.3')
%}
