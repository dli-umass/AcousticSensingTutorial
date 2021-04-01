%% System Configuration
total_num_of_chirps = 100;          % Total number of chirps needed to be processed
display_flag = 1;                   % Flag for displaying the result
display_flag_gt = 1;                % Flag for displaying the ground truth

paras.system_config.total_num_of_chirps = total_num_of_chirps;
paras.system_config.display_flag = display_flag;
paras.system_config.display_flag_gt = display_flag_gt;
%% FMCW Configuration
Vs = 343;                       % Sound speed
Fs = 44100;                     % Sampling frequency
Fc = 18000;                     % Start frequency
B = 4000;                       % Bandwidth
T = 0.04;                       % Sweep time
single_chirp_len = Fs*T;    	% Length of a single chirp

paras.fmcw_config.Vs = Vs;
paras.fmcw_config.Fs = Fs;
paras.fmcw_config.Fc = Fc;
paras.fmcw_config.B = B;
paras.fmcw_config.T = T; 
paras.fmcw_config.single_chirp_len = single_chirp_len;
paras.fmcw_config.total_time = 4;  % For simulation
%% Hardware Configuration
mic_int = 0.01070;                              % Microphone spacing
num_of_mics = 4;                
mic_int_vec = (0:(num_of_mics-1))*mic_int;      % Uniform array

paras.hardware_config.mic_int = mic_int;
paras.hardware_config.num_of_mics = num_of_mics;
paras.hardware_config.mic_int_vec = mic_int_vec;
%% Algorithm Configuration
%% ------ common ------
subsampling_factor = 40;      % set to 1 to avoid ,subsampling       
subsampling_indices = 1:subsampling_factor:single_chirp_len;
subsamp_single_chirp_len = ceil(single_chirp_len/subsampling_factor);
paras.algo_config.subsamp_factor = subsampling_factor;
paras.algo_config.subsamp_indices = subsampling_indices;
paras.algo_config.subsamp_single_chirp_len = subsamp_single_chirp_len ;
% ----- number of parameters ------
num_of_paras = 4;
paras.algo_config.num_of_paras = num_of_paras;
% ------ aoa ------ 
aoa_search_scope_min = 60;
aoa_search_scope_max = 140;
aoa_search_scope = [aoa_search_scope_min aoa_search_scope_max 2 20];                   % min; max; search stepsize; display stepsize
aoa_search = aoa_search_scope(1):aoa_search_scope(3):aoa_search_scope(2);
paras.algo_config.all_aoa_search_scope = aoa_search_scope;
paras.algo_config.all_aoa_search = aoa_search;
paras.algo_config.aoa_search_scope = aoa_search_scope;
paras.algo_config.aoa_search = aoa_search;
% ------ distance ------
dist_search_scope_min = 0;
dist_search_scope_max = 1;
dist_search_scope = [dist_search_scope_min dist_search_scope_max 0.01 0.2];            % min; max; search stepsize; display stepsize
dist_search = dist_search_scope(1):dist_search_scope(3):dist_search_scope(2);
paras.algo_config.all_dist_search_scope = dist_search_scope;
paras.algo_config.all_dist_search = dist_search;
paras.algo_config.dist_search_scope = dist_search_scope;
paras.algo_config.dist_search = dist_search;
% ------ velocity ------
num_of_chirps = 4;   % 2
vel_search_scope_min = -0.60;
vel_search_scope_max = 0.60;
vel_search_scope = [vel_search_scope_min vel_search_scope_max 0.05 0.02];
vel_search = vel_search_scope(1):vel_search_scope(3):vel_search_scope(2);
vel_ambiguity_interval = Vs/(2*(Fc+Fc+B)/2*T);
paras.algo_config.num_of_chirps = num_of_chirps;
paras.algo_config.all_vel_search = vel_search; 
paras.algo_config.all_vel_search_scope = vel_search_scope;
paras.algo_config.vel_search = vel_search;
paras.algo_config.vel_search_scope = vel_search_scope;
paras.algo_config.vel_ambiguity_interval = vel_ambiguity_interval;
%% ------ traditional estimation (FFT) ------
paras.algo_config.trad.dist_fft_size = 50*single_chirp_len;
paras.algo_config.trad.vel_fft_size = 50*num_of_chirps;
paras.algo_config.trad.aoa_fft_size = 100*num_of_mics;
%% ------ super-resolution estimation (MUSIC) ------
SUB_N = floor((single_chirp_len-subsampling_factor)/subsampling_factor)+1;              % number of samples after subsampling in each chirp
paras.algo_config.super.smoothed_window_dist = SUB_N-10;                                % set to SUB_N to avoid data smoothing
paras.algo_config.super.subsampling_factor_dist = subsampling_factor;
paras.algo_config.super.smoothed_window_vel = num_of_chirps;
paras.algo_config.super.smoothed_window_aoa = 4;