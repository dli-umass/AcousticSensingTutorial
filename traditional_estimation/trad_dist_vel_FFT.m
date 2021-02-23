function all_sig_path_paras = trad_dist_vel_FFT(sig,paras)
%TRAD_RANGE_VEL_FFT Estimate the range and velocity using FFT
%   sig             input signal
%   paras           parameters

%% Parameters
peak_thresh = 0.5;      % peak threshold: the amplitude of the peak must be larger than 0.5*max_peak_amp

Vs = paras.fmcw.Vs;
Fs = paras.fmcw.Fs;
Fc = paras.fmcw.Fc;
B = paras.fmcw.B;
T = paras.fmcw.T;

gt_mp_dists = paras.multipath.dists;
gt_mp_vels = paras.multipath.vels;
num_of_mps = paras.multipath.num_of_mps;

num_of_chirps = paras.algos.num_of_chirps;
display_flag = paras.algos.display_flag;
dist_max = paras.algos.dist_max;
vel_max = paras.algos.vel_max;
dist_fft_size = paras.algos.trad.dist_fft_size;
vel_fft_size = paras.algos.trad.vel_fft_size;
%% Algos
dist_search = linspace(0,Fs/2,dist_fft_size/2)*Vs*T/(2*B);
dist_idx = dist_search <= dist_max;
dist_search = dist_search(dist_idx);

doppler_freq = 1/T;
vel_search = -linspace(-doppler_freq/2,doppler_freq/2,vel_fft_size)*Vs/(2*Fc);
vel_idx = vel_search <=  vel_max;
vel_search = vel_search(vel_idx);

% ------ Range FFT ------
dist_chirp_fft = fft(sig,dist_fft_size,1);                          % perform column-wise fft
dist_chirp_fft = dist_chirp_fft(dist_idx,:);                        % choose the data within a certain distance, e.g. <= dist_max
% ------ Doppler FFT ------
dist_vel_fft = fftshift(fft(dist_chirp_fft,vel_fft_size,2),2);      % perform row-wise fft
dist_vel_fft = dist_vel_fft(:,vel_idx);                             % choose the data within a certain velocity, e.g. <= vel_max

dist_vel_map = abs(dist_vel_fft).^2/dist_fft_size/vel_fft_size;

% Estimate the range and velocity from the vel_dist_map
reg_max = imregionalmax(dist_vel_map);
[peak_row,peak_col,~] = find(reg_max);
peak_dist_vec = dist_search(peak_row);
peak_vel_vec = vel_search(peak_col);
peak_amp_vec = zeros(1, length(peak_row));
for peak_idx = 1:length(peak_row)
    peak_amp_vec(peak_idx) = dist_vel_map(peak_row(peak_idx),peak_col(peak_idx));
end
[~,peak_amp_idx] = sort(peak_amp_vec,'descend');
max_peak_amp = peak_amp_vec(peak_amp_idx(1));
all_sig_path_paras = cell(num_of_mps,1);
for mp_idx=1:num_of_mps
    if peak_amp_vec(peak_amp_idx(mp_idx)) >= peak_thresh*max_peak_amp
        sig_para_vec.raw_dist = peak_dist_vec(peak_amp_idx(mp_idx));
        sig_para_vec.vel = peak_vel_vec(peak_amp_idx(mp_idx));
        sig_para_vec.aoa = 0;
        sig_para_vec.amp = 0;
        % distance offse including range doppler coupling effect and movement
        dist_offset = -sig_para_vec.vel*Fc*T/B-sig_para_vec.vel*T*num_of_chirps/2;
        sig_para_vec.revised_dist = sig_para_vec.raw_dist - dist_offset;
        all_sig_path_paras{mp_idx} = sig_para_vec;
    else
        sig_para_vec.raw_dist = all_sig_path_paras{1}.raw_dist;
        sig_para_vec.vel = all_sig_path_paras{1}.vel;
        sig_para_vec.aoa = all_sig_path_paras{1}.aoa;
        sig_para_vec.amp = all_sig_path_paras{1}.amp;
        sig_para_vec.revised_dist = all_sig_path_paras{1}.revised_dist;
        all_sig_path_paras{mp_idx} = sig_para_vec;
    end
end

%% Display
if display_flag
    celldisp(all_sig_path_paras);
    
    % ------ 2D ------
    figure;
    ax = pcolor(dist_search,vel_search,dist_vel_map.');
%     points = zeros(2*num_of_mps,1);
%     for mp_idx=1:num_of_mps
%         hold on;
%         points(2*(mp_idx-1)+1) = plot(gt_mp_dists(mp_idx),gt_mp_vels(mp_idx),'+g','markersize',7,'markerfacecolor','g','linewidth',2);
%         hold on;
%         points(2*mp_idx) = plot(all_sig_path_paras{mp_idx}.revised_dist,all_sig_path_paras{mp_idx}.vel,'+r','markersize',7,'markerfacecolor','w','linewidth',2);
%     end
    points = zeros(num_of_mps,1);
    for mp_idx=1:num_of_mps
        hold on;
        points(mp_idx) = plot(gt_mp_dists(mp_idx),gt_mp_vels(mp_idx),'+g','markersize',7,'markerfacecolor','g','linewidth',2);
%         hold on;
%         points(2*mp_idx) = plot(all_sig_path_paras{mp_idx}.revised_dist,all_sig_path_paras{mp_idx}.vel,'+r','markersize',7,'markerfacecolor','w','linewidth',2);
    end
    set(ax, 'LineStyle','none');
    xlabel('Range (m)');
    xlim([0 dist_max]);
    xticks(0:0.2:dist_max);
    ylabel('Velocity(m/s)');
    yticks(-vel_max:0.02:vel_max);
    title('Range-Velocity FFT');
%     legend(points,'Groundtruth','Revised');
    legend(points,'Groundtruth');
    set(gca,'linewidth',1.5,'fontsize',20,'fontname','Times');
    colormap(jet);
    shading interp;
    colorbar;
    
    % ------ 3D ------
    figure;
    ax = surf(dist_search,vel_search,dist_vel_map.');
    for mp_idx=1:num_of_mps
        hold on;
        scatter3(gt_mp_dists(mp_idx),gt_mp_vels(mp_idx),peak_amp_vec(peak_amp_idx(mp_idx)),'MarkerEdgeColor','g','MarkerFaceColor','g');
        hold on;
        scatter3(all_sig_path_paras{mp_idx}.revised_dist,all_sig_path_paras{mp_idx}.vel,peak_amp_vec(peak_amp_idx(mp_idx)),'MarkerEdgeColor','r','MarkerFaceColor','r');
    end
    set(ax, 'LineStyle','none');
    xlabel('Range (m)');
    xlim([0 dist_max]);
    xticks(0:0.2:dist_max);
    ylabel('Velocity(m/s)');
    yticks(-vel_max:0.02:vel_max);
    zlabel('Amplitude');
    title('Range-Velocity FFT');
    set(gca,'linewidth',1.5,'fontsize',20,'fontname','Times');
    colormap(jet);
    shading interp;
    colorbar;
end

end

