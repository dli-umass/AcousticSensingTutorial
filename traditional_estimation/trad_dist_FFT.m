function all_sig_path_paras = trad_dist_FFT(sig,paras)
%TRAD_RANGE_FFT  Estimate the range using FFT
%   sig             input signal
%   paras           parameters

%% Parameters
Vs = paras.fmcw.Vs;
B = paras.fmcw.B;
T = paras.fmcw.T;
Fs = paras.fmcw.Fs;

dist_fft_size = paras.algos.trad.dist_fft_size;
dist_max = paras.algos.dist_max;
display_flag = paras.algos.display_flag;

num_of_mps = paras.multipath.num_of_mps;
gt_mp_dists = paras.multipath.dists;
%% FFT
dist_search = linspace(0,Fs/2,dist_fft_size/2)*Vs*T/(2*B);
dist_idx = dist_search <= dist_max;
dist_search = dist_search(dist_idx);

% ------ Range FFT ------
dist_fft = fft(sig,dist_fft_size,1);        % perform row-wise fft
dist_fft = dist_fft(dist_idx,:);            % choose the data within a certain distance, e.g. <= dist_max

dist_map = abs(dist_fft).^2/dist_fft_size;

chirp_idx = 1;
selected_chirp = dist_map(:,chirp_idx); 
[pks_amps, pks_locs] = findpeaks(selected_chirp);
[~, pks_idx] = sort(pks_amps,'descend');

max_peak_amp = pks_amps(pks_idx(1));
all_sig_path_paras = cell(1,num_of_mps);
sig_para_vec.dist = dist_search(pks_locs(pks_idx(1)));
sig_para_vec.amp = 0;
sig_para_vec.aoa = 0;
sig_para_vec.vel = 0;
all_sig_path_paras{1} = sig_para_vec;
for mp_idx=2:num_of_mps
    if pks_amps(pks_idx(mp_idx)) >= 0.5*max_peak_amp
        % remove the false positive samples (distance between two targets is a mulitplier of 0.5cm)
        if abs(gt_mp_dists(1)-all_sig_path_paras{1}.dist) >= 0.008
            sig_para_vec.dist = all_sig_path_paras{1}.dist;
        else
            sig_para_vec.dist =  dist_search(pks_locs(pks_idx(mp_idx)));
        end
    else
        sig_para_vec.dist = all_sig_path_paras{1}.dist;
    end
    sig_para_vec.amp = 0;
    sig_para_vec.aoa = 0;
    sig_para_vec.vel = 0;
    all_sig_path_paras{mp_idx} = sig_para_vec;
end
%% Display
if display_flag
    celldisp(all_sig_path_paras);
    
    figure;
    plot(dist_search,selected_chirp,'linewidth',3);
    hold on;
    point = plot([gt_mp_dists(1) gt_mp_dists(1)],[min(selected_chirp) max(selected_chirp)],'--r','linewidth',3);
    for mp_idx=2:num_of_mps
        hold on;
        plot([gt_mp_dists(mp_idx) gt_mp_dists(mp_idx)],[min(selected_chirp) max(selected_chirp)],'--r','linewidth',3);
    end
    xlabel('Range (m)');
    ylabel('Amplitude');
    title('Range FFT');
    legend(point,'Groundtruth');
    set(gca,'linewidth',1.5,'fontsize',20,'fontname','Times');
    
     % ------ 2D ------
%     figure;
%     ax = pcolor(dist_search,1:num_of_chirps,dist_map.');
%     for mp_idx=1:num_of_mps
%         hold on;
%         plot([gt_mp_dists(mp_idx) gt_mp_dists(mp_idx)],[0 num_of_chirps],'--w','linewidth',3);
%     end
%     set(ax, 'LineStyle','none');
%     xlabel('Range (m)');
%     xlim([0 dist_max]);
%     xticks(0:0.2:dist_max);
%     ylabel('Chirp Number');
%     title('Range-Chirp FFT');
%     set(gca,'linewidth',1.5,'fontsize',20,'fontname','Times');
%     colormap(jet);
%     shading interp;
%     colorbar;
    
    % ------ 3D ------
%     figure;
%     ax = surf(dist_search,1:num_of_chirps,dist_map.');
%     set(ax, 'LineStyle','none');
%     xlabel('Range (m)');
%     xlim([0 dist_max]);
%     xticks(0:0.2:dist_max);
%     ylabel('Chirp Number');
%     zlabel('Amplitude');
%     title('Distance-Chirp FFT');    
%     set(gca,'linewidth',1.5,'fontsize',20,'fontname','Times');
%     colormap(jet);
%     shading interp;
%     colorbar;
end

end

