function array_mix_sw = generate_mixed_sw(paras)
%GENERATEMIXEDSW Generate mixed soundwave
%   paras   parameters 

Vs = paras.fmcw.Vs;
Fs = paras.fmcw.Fs;
Fc = paras.fmcw.Fc;
B = paras.fmcw.B;
T = paras.fmcw.T;
single_chirp_len = paras.fmcw.single_chirp_len;

num_of_chirps = paras.algos.num_of_chirps;
num_of_mics = paras.algos.num_of_mics;
mic_int = paras.algos.mic_int;

amps = paras.multipath.amps;
aoas = paras.multipath.aoas;
dists = paras.multipath.dists;
vels = paras.multipath.vels;
num_of_mps = paras.multipath.num_of_mps;

%% Generate transmit and receive soundwave
amp_tx = 1; % amplitude for transmit soundwave 
init_phase = 0; % initial phase 
% ------ transmit soundwave ------
[trans_sw_cos,trans_sw_sin,~] = generate_transmit_sw(amp_tx,init_phase,Fs,Fc,B,T,T);
% ------ receive soundwave ------
array_rece_sw = zeros(num_of_mics,single_chirp_len*num_of_chirps);
for mic_idx=1:num_of_mics
    for mp_idx=1:num_of_mps
        prop_delay = (mic_idx-1)*mic_int*cosd(aoas(mp_idx))/Vs + 2*dists(mp_idx)/Vs;
        [mp_rece_sw,~] = generate_receive_sw(amps(mp_idx),init_phase,vels(mp_idx),prop_delay,paras);
        array_rece_sw(mic_idx, :) = squeeze(array_rece_sw(mic_idx, :)) + mp_rece_sw(1:single_chirp_len*num_of_chirps);
    end
    % ------ add gaussian noise ------
%     array_rece_sw(mic_idx, :) = awgn(array_rece_sw(mic_idx, :),3,'measured');
end
%% Construct mixed soundwave
array_mix_sw = zeros(single_chirp_len,num_of_mics,num_of_chirps);
% filter out the received soundwave whose frequency is between Fc and Fc +B
bpFilt_chirp = fir1(80,[(Fc-1000)/(Fs/2) (Fc+B+1000)/(Fs/2)],'bandpass');
% filter out the low frequency part
lpFilt_chirp = fir1(80,[(1)/(Fs/2) (1000)/(Fs/2)],'bandpass');
for chirp_idx=1:num_of_chirps
    for mic_idx=1:num_of_mics
        rece_sw = array_rece_sw(mic_idx,(chirp_idx-1)*single_chirp_len+1:chirp_idx*single_chirp_len);
        rece_sw = filtfilt(bpFilt_chirp,1,rece_sw);
        
        % signal pluralizing
        mix_sw_cos = rece_sw.*trans_sw_cos;
        mix_sw_sin = rece_sw.*trans_sw_sin;
        
        mix_sw_cos = filtfilt(lpFilt_chirp,1,mix_sw_cos);
        mix_sw_sin = filtfilt(lpFilt_chirp,1,mix_sw_sin);
        
        array_mix_sw(:,mic_idx,chirp_idx) = mix_sw_cos + 1j*mix_sw_sin;
    end
end

end

