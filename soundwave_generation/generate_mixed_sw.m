function array_mix_sw = generate_mixed_sw(paras)
%GENERATEMIXEDSW Generate mixed soundwave
%   paras   parameters 
%% Parameters
Vs = paras.fmcw_config.Vs;
Fs = paras.fmcw_config.Fs;
Fc = paras.fmcw_config.Fc;
B = paras.fmcw_config.B;
single_chirp_len = paras.fmcw_config.single_chirp_len;

total_num_of_chirps = paras.system_config.total_num_of_chirps;

num_of_mics = paras.hardware_config.num_of_mics;
mic_int_vec = paras.hardware_config.mic_int_vec;

amps = paras.multipath.gt_amps;
aoas = paras.multipath.gt_aoas;
dists = paras.multipath.gt_dists;
vels = paras.multipath.gt_vels;
num_of_mps = paras.multipath.num_of_mps;
%% Generate transmit and receive soundwave
amp_tx = 10;         % amplitude for transmit soundwave 
init_phase = 0;      % initial phase 
% ------ transmit soundwave ------
[trans_sw_cos,trans_sw_sin,~] = generate_transmit_sw(amp_tx,init_phase,paras);
% ------ middle soundwave ------
% [middle_sw_cos,middle_sw_sin,~] = generate_middle_sw(amp_tx,init_phase,paras);
% ------ receive soundwave ------
array_rece_sw = zeros(num_of_mics,single_chirp_len*total_num_of_chirps);
for mic_idx=1:num_of_mics
    for mp_idx=1:num_of_mps
        prop_delay = mic_int_vec(mic_idx)*cosd(aoas(mp_idx))/Vs + 2*dists(mp_idx)/Vs;
        [mp_rece_sw,~] = generate_receive_sw(amps(mp_idx),init_phase,vels(mp_idx),prop_delay,paras);
%         [mp_rece_sw,rece_freq] = generate_receive_sw(amps(mp_idx),init_phase,vels(mp_idx),prop_delay,paras);
%         figure;
%         plot(rece_freq);
        array_rece_sw(mic_idx,:) = squeeze(array_rece_sw(mic_idx,:)) + mp_rece_sw(1:single_chirp_len*total_num_of_chirps);
    end
    % ------ add gaussian noise ------
%     array_rece_sw(mic_idx, :) = awgn(array_rece_sw(mic_idx, :),0,'measured');
end
%% Construct mixed soundwave
array_mix_sw = zeros(single_chirp_len,num_of_mics,total_num_of_chirps);
% filter out the received soundwave whose frequency is between Fc and Fc +B
bpFilt_chirp = fir1(80,[(Fc-10)/(Fs/2) (Fc+B+10)/(Fs/2)],'bandpass');
% filter out the low frequency part
lpFilt_chirp = fir1(80,[(1)/(Fs/2) (1000)/(Fs/2)],'bandpass');
for chirp_idx=1:total_num_of_chirps
    for mic_idx=1:num_of_mics
       %% Apply a band-pass filter to remove out-of-band noise
        rece_sw = array_rece_sw(mic_idx,(chirp_idx-1)*single_chirp_len+1:chirp_idx*single_chirp_len).';
        rece_sw = filtfilt(bpFilt_chirp,1,rece_sw);
       %% Multiply the middle and received chirps
        mix_sw_cos = rece_sw.*trans_sw_cos;
        mix_sw_sin = rece_sw.*trans_sw_sin;
       %% Apply a low-pass filter to obtain desired sinusoids
        mix_sw_cos = filtfilt(lpFilt_chirp,1,mix_sw_cos);
        mix_sw_sin = filtfilt(lpFilt_chirp,1,mix_sw_sin);
       %% Signal Pluralizing
        array_mix_sw(:,mic_idx,chirp_idx) = mix_sw_cos + 1j*mix_sw_sin;
    end
end
end

