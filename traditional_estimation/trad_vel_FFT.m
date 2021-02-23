function all_sig_path_paras = trad_vel_FFT(sig,paras)
%TRAD_VEL_FFT  Estimate the velocity using FFT
%   sig             input signal
%   paras           parameters

%% Parameters
Vs = paras.fmcw.Vs;
T = paras.fmcw.T;
Fc = paras.fmcw.Fc;
single_chirp_len = paras.fmcw.single_chirp_len;

vel_fft_size = paras.algos.trad.vel_fft_size;
vel_max = paras.algos.vel_max;
display_flag = paras.algos.display_flag;

num_of_mps = paras.multipath.num_of_mps;
gt_mp_vels = paras.multipath.vels;
%% Algos
doppler_freq = 1/T;
vel_search = -linspace(-doppler_freq/2,doppler_freq/2,vel_fft_size)*Vs/(2*Fc);
vel_idx = vel_search <=  vel_max;
vel_search = vel_search(vel_idx);

% ------ Doppler FFT ------
vel_fft = fftshift(fft(sig,vel_fft_size,2),2);      % perform column-wise fft
vel_fft = vel_fft(:,vel_idx);                       % choose the data within a certain velocity, e.g. <= vel_max

vel_map = abs(vel_fft).^2/vel_fft_size;

sample_idx = 100;
selected_sample = vel_map(sample_idx,:);
[pks_amps,pks_locs] = findpeaks(selected_sample);
[~, pks_idx] = sort(pks_amps,'descend');

max_peak_amp = pks_amps(pks_idx(1));
all_sig_path_paras = cell(1,num_of_mps);
sig_para_vec.dist = 0;
sig_para_vec.amp = 0;
sig_para_vec.aoa = 0;
sig_para_vec.vel = vel_search(pks_locs(pks_idx(1)));
all_sig_path_paras{1} = sig_para_vec;
for mp_idx=2:num_of_mps
    if pks_amps(pks_idx(mp_idx)) >= 0.5*max_peak_amp
        sig_para_vec.vel =  vel_search(pks_locs(pks_idx(mp_idx)));
    else
        sig_para_vec.vel = all_sig_path_paras{1}.vel;
    end
    sig_para_vec.dist = 0;
    sig_para_vec.amp = 0;
    sig_para_vec.aoa = 0;
    all_sig_path_paras{mp_idx} = sig_para_vec;
end
%% Display
if display_flag
    celldisp(all_sig_path_paras);
    
    % ------ 1D ------
    figure;
    plot(vel_search,selected_sample,'linewidth',3);
    for mp_idx=1:num_of_mps
        hold on;
        plot([gt_mp_vels(mp_idx) gt_mp_vels(mp_idx)],[min(selected_sample) max(selected_sample)],'--r','linewidth',3);
    end
    xlabel('Velocity (m/s)');
    ylabel('Amplitude');
    title('Velocity FFT');
    set(gca,'linewidth',1.5,'fontsize',20,'fontname','Times');
    
    % ------ 2D ------
    figure;
    ax = pcolor(1:single_chirp_len,vel_search,vel_map.');
    for mp_idx=1:num_of_mps
        hold on;
        plot([1 single_chirp_len],[gt_mp_vels(mp_idx) gt_mp_vels(mp_idx)],'--w','linewidth',3);
    end
    set(ax, 'LineStyle','none');
    xlabel('Sample Number');
    ylabel('Velocity(m/s)');
    yticks(-vel_max:0.02:vel_max);
    title('Samples-Velocity FFT');
    set(gca,'linewidth',1.5,'fontsize',20,'fontname','Times');
    colormap(jet);
    shading interp;
    colorbar;
    
    % ------ 3D ------
%     figure;
%     ax = surf(1:single_chirp_len,vel_search,vel_map.');
%     set(ax, 'LineStyle','none');
%     xlabel('Sample Number');
%     ylabel('Velocity(m/s)');
%     yticks(-vel_max:0.02:vel_max);
%     zlabel('Amplitude');
%     title('Samples-Velocity FFT');
%     set(gca,'linewidth',1.5,'fontsize',20,'fontname','Times');
%     colormap(jet);
%     shading interp;
%     colorbar;
end

end

