clc;
clear;
close all;
%% Add path
addpath('soundwave_generation');
addpath('traditional_estimation');
addpath('super_estimation');
%% Global Configuration
global_config;
%% Parameters for the multipath
% ------ One Target ------
% aoas = 100;
% dists = 0.44;
% amps = 0.0008/dists(1)^2;
% vels = 0;

% ------ Two Targets ------
aoas = [83 100];
dists = [0.41 0.44];
amps = [0.001/dists(1)^2 0.0008/dists(2)^2]; 
vels = [0 0];

num_of_mps = length(amps);
paras.multipath.gt_amps = amps;
paras.multipath.gt_aoas = aoas;
paras.multipath.gt_dists = dists;
paras.multipath.gt_vels = vels;
paras.multipath.revised_dists = dists;
paras.multipath.num_of_mps = num_of_mps;
%% Generate Mixed Soundwave
array_mix_sw = generate_mixed_sw(paras);
%% Traditional Estimation (FFT)
%% ------ distance estimation -----
fprintf('------ distance estimation by FFT ------\n');
mix_sw = squeeze(array_mix_sw(:,1,1:num_of_chirps));
trad_dist_FFT(mix_sw, paras);
%% ------ distance + aoa estimation -----
fprintf('------ distance + AoA estimation by FFT ------\n');
mix_sw = squeeze(array_mix_sw(:,:,1));
trad_dist_aoa_FFT(mix_sw, paras);
%% Super-resolution Estimation
%% ------ distance estimation -----
fprintf('------ distance estimation by MUSIC ------\n');
mix_sw = squeeze(array_mix_sw(:,1,1));
super_dist(mix_sw, paras);
%% ------ distance + aoa estimation -----
fprintf('------ distance + AoA estimation by 2D MUSIC ------\n');
mix_sw = squeeze(array_mix_sw(:,:,1));
super_dist_aoa(mix_sw.', paras);