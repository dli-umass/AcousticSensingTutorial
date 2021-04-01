clc;
clear;
close all;
%% Add path
addpath('soundwave_generation');
addpath('traditional_estimation');
addpath('super_estimation');
%% Parametes for FMCW signals
Vs = 343;                   % Sound velocity
Fs = 44100;                 % Sampling frequency
Fc = 16000;                 % Start frequency
B = 4000;                   % Bandwidth
T = 0.04;                    % Sweep time
total_time = 2;             % Total time
single_chirp_len = Fs*T;    % Length of a single chirp

paras.fmcw.Vs = Vs;
paras.fmcw.Fs = Fs;
paras.fmcw.Fc = Fc;
paras.fmcw.B = B;
paras.fmcw.T = T;
paras.fmcw.single_chirp_len = single_chirp_len;
paras.fmcw.total_time = total_time;
%% Parameters for the multipath
aoas = [90 80];
dists = [0.3 0.31];
amps = [0.001/dists(1)^2 0.001/dists(2)^2]; 
vels = [0 0.01];

% aoas = [90];
% dists = [0.2];
% amps = [0.001/dists(1)^4];
% vels = [-0.04];

num_of_mps = length(amps);
paras.multipath.amps = amps;
paras.multipath.aoas = aoas;
paras.multipath.dists = dists;
paras.multipath.vels = vels;
paras.multipath.num_of_mps = num_of_mps;
%% Common parameters for algos
% ------ flag for display ------
display_flag = true;
paras.algos.display_flag = display_flag;
paras.algos.display_flag_gt = display_flag_gt; 

wave_len = Vs/(Fc+B);
% ------ aoa ------
num_of_mics = 4;
mic_int = 1/2*wave_len;
aoa_max = 180;
paras.algos.num_of_mics = num_of_mics;
paras.algos.mic_int = mic_int;
paras.algos.aoa_max = aoa_max;
% ------ distacne ------
dist_max = 1;
paras.algos.dist_max = dist_max;
% ------ velocity ------
num_of_chirps = 20;
vel_max = 0.1;
paras.algos.num_of_chirps = num_of_chirps;
paras.algos.vel_max = vel_max;
%% Parameters for traditional estimation
% ------ FFT -----
paras.algos.trad.dist_fft_size = 50*single_chirp_len;
paras.algos.trad.vel_fft_size = 10*num_of_chirps;
paras.algos.trad.aoa_fft_size = 50*num_of_mics;
%% Generate Mixed Soundwave
array_mix_sw = generate_mixed_sw(paras);
%% Traditional Estimation
fprintf('------ Traditional Estimation  ------\n');
% ------ Distance ------
fprintf('- FFT Distance -\n');
mix_sw = squeeze(array_mix_sw(:,1,:));
trad_dist_FFT(mix_sw,paras);
% ------ Distance + Velocity ------
fprintf('- FFT Distance -\n');
mix_sw = squeeze(array_mix_sw(:,4,:));
trad_dist_vel_FFT(mix_sw,paras);
%% Super-resolution Estimation
% ------ distance estimation -----
fprintf('------ distance estimation by MUSIC ------\n');
mix_sw = squeeze(array_mix_sw(:,1,1));
super_dist(mix_sw, paras);
%% ------ distance + aoa estimation -----
fprintf('------ distance + AoA estimation by 2D MUSIC ------\n');
mix_sw = squeeze(array_mix_sw(:,:,1));
super_dist_aoa(mix_sw.', paras);