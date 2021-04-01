function all_sig_path_paras = super_dist(sig,paras)
%TRAD_DIST_MUSIC Estimate the range using MUSIC
%   sig                     input signal
%   paras                   parameters
%   all_sig_path_paras      path parameters for all signals (cell)
%% Parameters
peak_thresh = 0;       % peak threshold: the amplitude of the peak must be larger than peak_thresh*max_peak_amp

display_flag = paras.system_config.display_flag;
display_flag_gt = paras.system_config.display_flag_gt;

Vs = paras.fmcw_config.Vs;
B = paras.fmcw_config.B;
T = paras.fmcw_config.T;
Fs = paras.fmcw_config.Fs;

dist_search_scope = paras.algo_config.dist_search_scope;
subsampling_factor = paras.algo_config.subsamp_factor;
smooth_window = paras.algo_config.super.smoothed_window_dist;

gt_mp_dists = paras.multipath.gt_dists;
num_of_mps = paras.multipath.num_of_mps;
%% Algorithm
dist_min = dist_search_scope(1);
dist_max = dist_search_scope(2);
dist_step = dist_search_scope(3);
freq_min = 2*dist_min*B/(T*Vs);         % minimum frequency 
freq_max = 2*dist_max*B/(T*Vs);         % maximum frequency
freq_step = 2*dist_step*B/(T*Vs);        % frequency stepsize
[fom_freq,fom_sp] = music_1D(sig,Fs,freq_min,freq_max,freq_step,subsampling_factor,smooth_window);
dist_search = fom_freq/B*T*Vs/2;

fom_sp = abs(fom_sp);
fom_sp = fom_sp/max(fom_sp);
[pks_amps, pks_locs] = findpeaks(fom_sp);
[~, pks_idx] = sort(pks_amps,'descend');

max_peak_amp = pks_amps(pks_idx(1));
all_sig_path_paras = cell(1,num_of_mps);
for mp_idx=1:num_of_mps
    if mp_idx <= length(pks_idx)
        if pks_amps(pks_idx(mp_idx)) >= peak_thresh*max_peak_amp
            est_dist = dist_search(pks_locs(pks_idx(mp_idx)));
        else
            % avoid side lobes
            est_dist = all_sig_path_paras{1}.raw_dist;     
        end
    else 
        % avoid the case where # estimated paths are less than # true paths
        non_empty_idx = find(~cellfun(@isempty, all_sig_path_paras));
        last_idx = non_empty_idx(end);
        est_dist = all_sig_path_paras{last_idx}.raw_dist;
    end

    sig_para_vec.amp = 0;
    sig_para_vec.raw_dist = est_dist;
    sig_para_vec.raw_aoa = 0;
    sig_para_vec.raw_vel = 0;
    all_sig_path_paras{mp_idx} = sig_para_vec;
end
%% Display
if display_flag
    celldisp(all_sig_path_paras);
    
    dist_display_step = dist_search_scope(4);
    figure;
    plot(dist_search,fom_sp.','linewidth',3);
    if display_flag_gt
        % ------ Ground Truth -----
        hold on;
        point = plot([gt_mp_dists(1) gt_mp_dists(1)],[min(pks_amps) max(pks_amps)],'--r','linewidth',3);
        for mp_idx=2:num_of_mps
            hold on;
            plot([gt_mp_dists(mp_idx) gt_mp_dists(mp_idx)],[min(pks_amps) max(pks_amps)],'--r','linewidth',3);
        end
        legend(point,'Groundtruth');
    end
    xlabel('Range (m)');
    xlim([dist_min dist_max]);
    xticks(dist_min:dist_display_step:dist_max);
    ylabel('Amplitude');
    title('Range MUSIC');
    set(gca,'linewidth',1.5,'fontsize',20,'fontname','Arial');
end
end

