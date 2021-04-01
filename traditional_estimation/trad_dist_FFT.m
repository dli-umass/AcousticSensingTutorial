function [all_sig_path_paras,selected_chirp] = trad_dist_FFT(sig,paras)
%TRAD_RANGE_FFT  Estimate the range using FFT
%   sig             input signal
%   paras           parameters

%% Parameters
peak_thresh = 0;       % peak threshold: the amplitude of the peak must be larger than peak_thresh*max_peak_amp

Vs = paras.fmcw_config.Vs;
Fs = paras.fmcw_config.Fs;
B = paras.fmcw_config.B;
T = paras.fmcw_config.T;

display_flag = paras.system_config.display_flag;
display_flag_gt = paras.system_config.display_flag_gt;

dist_search_scope = paras.algo_config.dist_search_scope;
num_of_chirps = paras.algo_config.num_of_chirps;
dist_fft_size = paras.algo_config.trad.dist_fft_size;

num_of_mps = paras.multipath.num_of_mps;
gt_mp_dists = paras.multipath.gt_dists;
%% Algorithm
dist_min = dist_search_scope(1);
dist_max = dist_search_scope(2);
dist_search = linspace(0,Fs/2,dist_fft_size/2)*Vs*T/(2*B);
dist_idx = (dist_search >= dist_min) & (dist_search <= dist_max);
dist_search = dist_search(dist_idx);

% ------ Range FFT ------
dist_fft = fft(sig,dist_fft_size,1);        % perform row-wise fft
% dist_fft_2 = dist_fft(1:dist_fft_size/2+1);
% dist_fft_2(2:end-1) = flip(dist_fft(dist_fft_size/2+2:end));
% dist_fft(2:dist_fft_size/2,:) = flip(dist_fft(dist_fft_size/2+2:end,:),1);
dist_fft = dist_fft(dist_idx,:);            % choose the data within a certain distance, e.g. <= dist_max

dist_map = abs(dist_fft).^2/dist_fft_size;

chirp_idx = 1;
selected_chirp = dist_map(:,chirp_idx); 
% ------ Background subtraction test ------
% selected_chirp = dist_map(:,chirp_idx) - paras.back_fft; 
[pks_amps, pks_locs] = findpeaks(selected_chirp);
[~, pks_idx] = sort(pks_amps,'descend');

max_peak_amp = pks_amps(pks_idx(1));
all_sig_path_paras = cell(1,num_of_mps);
for mp_idx=1:num_of_mps
    if mp_idx <= length(pks_idx)
        if pks_amps(pks_idx(mp_idx)) >= peak_thresh*max_peak_amp
            est_dist =  dist_search(pks_locs(pks_idx(mp_idx)));
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
    plot(dist_search,selected_chirp,'linewidth',3);
    if display_flag_gt
        % ------ Ground Truth -----
        hold on;
        point = plot([gt_mp_dists(1) gt_mp_dists(1)],[min(selected_chirp) max(selected_chirp)],'--r','linewidth',3);
        for mp_idx=2:num_of_mps
            hold on;
            plot([gt_mp_dists(mp_idx) gt_mp_dists(mp_idx)],[min(selected_chirp) max(selected_chirp)],'--r','linewidth',3);
        end
        legend(point,'Groundtruth');
    end
    xlabel('Range (m)');
    xlim([dist_min dist_max]);
    xticks(dist_min:dist_display_step:dist_max);
    ylabel('Amplitude');
    title('Range FFT');
    set(gca,'linewidth',1.5,'fontsize',20,'fontname','Arial');
    
%      % ------ 2D ------
%     figure;
%     ax = pcolor(dist_search,1:num_of_chirps,dist_map.');
%     if display_flag_gt
%         % ------ Ground Truth -----
%         for mp_idx=1:num_of_mps
%             hold on;
%             plot([gt_mp_dists(mp_idx) gt_mp_dists(mp_idx)],[0 num_of_chirps],'--w','linewidth',3);
%         end
%     end
%     set(ax, 'LineStyle','none');
%     xlabel('Range (m)');
%     xlim([0 dist_max]);
%     xticks(0:dist_display_step:dist_max);
%     ylabel('Chirp Number');
%     title('Range-Chirp FFT');
%     set(gca,'linewidth',1.5,'fontsize',20,'fontname','Arial');
%     colormap(jet);
%     shading interp;
%     colorbar;
%     
%     % ------ 3D ------
%     figure;
%     ax = surf(dist_search,1:num_of_chirps,dist_map.');
%     set(ax, 'LineStyle','none');
%     xlabel('Range (m)');
%     xlim([0 dist_max]);
%     xticks(0:dist_display_step:dist_max);
%     ylabel('Chirp Number');
%     zlabel('Amplitude');
%     title('Distance-Chirp FFT');    
%     set(gca,'linewidth',1.5,'fontsize',20,'fontname','Arial');
%     colormap(jet);
%     shading interp;
%     colorbar;
end

end

