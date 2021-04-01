function [all_sig_path_paras] = trad_dist_aoa_FFT(sig,paras)
%TRAD_DIST_AOA_FFT Estimate the range and aoa using FFT
%   sig                     input signal
%   paras                   parameters
%   all_sig_path_paras      path parameters for all signals (cell)
%% Parameters
peak_thresh = 0;     % peak threshold: the amplitude of the peak must be larger than 0.5*max_peak_amp

display_flag = paras.system_config.display_flag;
display_flag_gt = paras.system_config.display_flag_gt;

Vs = paras.fmcw_config.Vs;
B = paras.fmcw_config.B;
T = paras.fmcw_config.T;
Fs = paras.fmcw_config.Fs;

% ------ distance ------
dist_search_scope = paras.algo_config.dist_search_scope;
dist_fft_size = paras.algo_config.trad.dist_fft_size;
% ------ aoa ------
aoa_search_scope = paras.algo_config.aoa_search_scope;
aoa_fft_size = paras.algo_config.trad.aoa_fft_size;

gt_mp_dists = paras.multipath.gt_dists;
gt_mp_aoas = paras.multipath.gt_aoas;
num_of_mps = paras.multipath.num_of_mps;
%% Algorithm
dist_min = dist_search_scope(1);
dist_max = dist_search_scope(2);
dist_search = linspace(0,Fs/2,dist_fft_size/2)*Vs*T/(2*B);
dist_idx = (dist_search >= dist_min) & (dist_search <= dist_max);
dist_search = dist_search(dist_idx);

aoa_min = aoa_search_scope(1);
aoa_max = aoa_search_scope(2);
aoa_search = acosd(linspace(-1,1,aoa_fft_size));
aoa_idx = (aoa_search>=aoa_min) & (aoa_search<=aoa_max);
aoa_search = aoa_search(aoa_idx);

% ------ Range FFT ------
dist_mic_fft = fft(sig,dist_fft_size,1);                            % perform column-wise fft
dist_mic_fft = dist_mic_fft(dist_idx,:);                            % choose the data within a certain distance
% ------ AoA FFT ------
dist_aoa_fft = fftshift(fft(dist_mic_fft,aoa_fft_size,2),2);        % perform row-wise fft
dist_aoa_fft = dist_aoa_fft(:,aoa_idx);                             % choose the data within a certain angle

dist_aoa_map = abs(dist_aoa_fft).^2/dist_fft_size/aoa_fft_size;

% ------ Estimate the range and velocity from the vel_dist_map ------
reg_max = imregionalmax(dist_aoa_map);
[peak_row,peak_col,~] = find(reg_max);
peak_dist_vec = dist_search(peak_row);
peak_aoa_vec = aoa_search(peak_col);
peak_amp_vec = zeros(1, length(peak_row));
for peak_idx = 1:length(peak_row)
    peak_amp_vec(peak_idx) = dist_aoa_map(peak_row(peak_idx),peak_col(peak_idx));
end
[~,peak_amp_idx] = sort(peak_amp_vec,'descend');
max_peak_amp = peak_amp_vec(peak_amp_idx(1));
all_sig_path_paras = cell(1,num_of_mps);
for mp_idx=1:num_of_mps
    if mp_idx <= length(peak_amp_idx)
        if peak_amp_vec(peak_amp_idx(mp_idx)) >= peak_thresh*max_peak_amp
            est_dist = peak_dist_vec(peak_amp_idx(mp_idx));
            est_aoa = peak_aoa_vec(peak_amp_idx(mp_idx));
        else
            % avoid side lobes
            est_dist = all_sig_path_paras{1}.raw_dist;
            est_aoa = all_sig_path_paras{1}.raw_aoa;
        end
    else
        % avoid the case where # estimated paths are less than # true paths
        est_dist = all_sig_path_paras{1}.raw_dist;
        est_aoa = all_sig_path_paras{1}.raw_aoa;
    end
    
    sig_para_vec.amp = 0;
    sig_para_vec.raw_dist = est_dist;
    sig_para_vec.raw_aoa = est_aoa;
    sig_para_vec.raw_vel = 0;
    all_sig_path_paras{mp_idx} = sig_para_vec;
end
%% Display
if display_flag
    celldisp(all_sig_path_paras);
    
    dist_display_step = dist_search_scope(4);
    aoa_display_step = aoa_search_scope(4);
    % ------ 2D ------
    figure;
    ax = pcolor(dist_search,aoa_search,dist_aoa_map.');
    if display_flag_gt
        % groundtruth
        points = zeros(1,1);
        hold on;
        points(1) = scatter(gt_mp_dists(1),gt_mp_aoas(1),80,'+r','linewidth',3);
        for mp_idx=2:num_of_mps
            hold on;
            scatter(gt_mp_dists(mp_idx),gt_mp_aoas(mp_idx),80,'+r','linewidth',3);
        end
        legend(points,'Groundtruth');
    end
    set(ax, 'LineStyle','none');
    xlabel('Range (m)');
    xlim([dist_min dist_max]);
    xticks(dist_min:dist_display_step:dist_max);
    ylabel('AoA (Deg)');
    ylim([aoa_min aoa_max]);
    yticks(aoa_min:aoa_display_step:aoa_max);
    title('Range-AoA FFT');
    set(gca,'linewidth',1.5,'fontsize',20,'fontname','Arial');
    shading interp;
    colorbar;
    
    % ------ 3D ------
%     figure;
%     ax = surf(dist_search,aoa_search,dist_aoa_map.');
%     if display_flag_gt
%         % groundtruth
%         points = zeros(1,1);
%         hold on;
%         points(1) = scatter3(gt_mp_dists(1),gt_mp_aoas(1),peak_amp_vec(peak_amp_idx(1)),'o','MarkerEdgeColor',[0.8500, 0.3250, 0.0980],'markerfacecolor','w','linewidth',2);
%         for mp_idx=2:num_of_mps
%             hold on;
%             scatter3(gt_mp_dists(mp_idx),gt_mp_aoas(mp_idx),peak_amp_vec(peak_amp_idx(mp_idx)),'o','MarkerEdgeColor',[0.8500, 0.3250, 0.0980],'markerfacecolor','w','linewidth',2);
%         end
%         legend(points,'Groundtruth','location','best');
%     end
%     set(ax, 'LineStyle','none');
%     xlabel('Range (m)');
%     xlim([dist_min dist_max]);
%     xticks(dist_min:dist_display_step:dist_max);
%     ylabel('AoA (Deg)');
%     ylim([aoa_min aoa_max]);
%     yticks(aoa_min:aoa_display_step:aoa_max);
%     zlabel('Amplitude');
%     title('Range-AoA FFT');
%     set(gca,'linewidth',1.5,'fontsize',20,'fontname','Arial');
%     colormap(jet);
%     shading interp;
%     colorbar;
end
end

