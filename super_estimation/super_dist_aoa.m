function all_sig_path_paras = super_dist_aoa(sig,paras)
%SUPER_DIST_AOA Estimate the distance and AoA simultaneously
%   sig             input signal
%   paras           parameters
%   Reference: "2D-MUSIC TECHNIQUE APPLIED TO A COHERENT FMCW MIMO RADAR"
%% Parameters
peak_thresh = 1e-6;      % peak threshold: the amplitude of the peak must be larger than 0.5*max_peak_amp

Vs = paras.fmcw_config.Vs;
Fs = paras.fmcw_config.Fs;
Fc = paras.fmcw_config.Fc;
B = paras.fmcw_config.B;
T = paras.fmcw_config.T;

display_flag = paras.system_config.display_flag;
display_flag_gt = paras.system_config.display_flag_gt;

gt_mp_dists = paras.multipath.gt_dists;
gt_mp_aoas = paras.multipath.gt_aoas;
num_of_mps = paras.multipath.num_of_mps;

mic_int_vec = paras.hardware_config.mic_int_vec;

% ------ distance ------
dist_search = paras.algo_config.dist_search;
dist_search_scope = paras.algo_config.dist_search_scope;
% ------ aoa -----
aoa_search = paras.algo_config.aoa_search;
aoa_search_scope = paras.algo_config.aoa_search_scope;

smoothed_window_aoa = paras.algo_config.super.smoothed_window_aoa;
smoothed_window_dist = paras.algo_config.super.smoothed_window_dist;
subsampling_factor_dist = paras.algo_config.super.subsampling_factor_dist;
%% Algorithm
Ts = 1/Fs;
[M,N] = size(sig);
sub_N = floor((N-subsampling_factor_dist)/subsampling_factor_dist)+1;
p1 = M-smoothed_window_aoa+1;
p2 = sub_N-smoothed_window_dist+1;

trans_matrix = eye(smoothed_window_aoa*smoothed_window_dist, smoothed_window_aoa*smoothed_window_dist);
J = fliplr(trans_matrix);

R_sub = 0;
X = zeros(smoothed_window_aoa*smoothed_window_dist, p1*p2);
if sub_N <= smoothed_window_dist
    fprintf('Warning: The subsampling factor K is too large!\n');
    return;
end

% ------ sub-sampling & data smoothing ------
for i=1:subsampling_factor_dist
    % ------ sub-sampling ------
    sub_sig = zeros(M,sub_N);
    sub_index = i:subsampling_factor_dist:i+subsampling_factor_dist*(sub_N-1);
    sub_sig(:,:) = sig(:,sub_index);
    
    % ------ smooth data ------
    for j = 1:p1
        for k = 1:p2
            X(:,(j-1)*p2+k) = reshape(sub_sig(j:j+smoothed_window_aoa-1,k:k+smoothed_window_dist-1).',[smoothed_window_aoa*smoothed_window_dist,1]);
        end
    end
    R_sub = R_sub + 1/(2*p1*p2)*(X*X'+J*(X*X').'*J);
end
R_sub = R_sub/subsampling_factor_dist;

[V,D] = eig(R_sub);                     % Compute the eigenvalues and eigenvectors of the auto-correlation matrix
[D_val,I] = sort(diag(D),'descend');    % Sort the eigenvectors in a descending order in terms of the magnitude of corresponding eigenvalues.

Level=10^4;
if abs(max(D_val)/min(D_val)) < Level
    fprintf('Warning: The eigen values may not correct!\n         The estinated DOAs may wrong!\n')
end
noise_num = smoothed_window_aoa*smoothed_window_dist-length(find(D_val<max(D_val)/Level));   % Compute the number of noises in the noise space

aoa_dist_map = zeros(length(aoa_search),length(dist_search));
time_vec = 0:subsampling_factor_dist*Ts:(smoothed_window_dist-1)*subsampling_factor_dist*Ts;
for aoa_idx=1:length(aoa_search)
    for dist_idx = 1:length(dist_search)
        tao = 2*dist_search(dist_idx)/Vs+mic_int_vec.'*cosd(aoa_search(aoa_idx))/Vs;
        steering_vector = exp(1i*2*pi*(Fc*tao+B/T*tao*time_vec));
        steering_vector = reshape(steering_vector.',[smoothed_window_aoa*smoothed_window_dist,1]);
        Rn = V(:, I(noise_num+1:smoothed_window_aoa*smoothed_window_dist));                                 % Noise space matrix
        aoa_dist_map(aoa_idx,dist_idx)=1/(steering_vector'*(Rn*Rn')*steering_vector);
    end
end
aoa_dist_map = abs(aoa_dist_map).^2;

% Estimate the range and aoa from the aoa_dist_map
reg_max = imregionalmax(aoa_dist_map);
[peak_row,peak_col,~] = find(reg_max);
peak_aoa_vec = aoa_search(peak_row);
peak_dist_vec = dist_search(peak_col);
peak_amp_vec = zeros(1, length(peak_row));
for peak_idx = 1:length(peak_row)
    peak_amp_vec(peak_idx) = aoa_dist_map(peak_row(peak_idx),peak_col(peak_idx));
end
[~,peak_amp_idx] = sort(peak_amp_vec,'descend');
max_peak_amp = peak_amp_vec(peak_amp_idx(1));
all_sig_path_paras = cell(num_of_mps,1);
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
        non_empty_idx = find(~cellfun(@isempty, all_sig_path_paras));
        last_idx = non_empty_idx(end);
        est_dist = all_sig_path_paras{last_idx}.raw_dist;
        est_aoa = all_sig_path_paras{last_idx}.raw_aoa;
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
    
    dist_min = dist_search_scope(1);
    dist_max = dist_search_scope(2);
    dist_display_step = dist_search_scope(4);
    aoa_min = aoa_search_scope(1);
    aoa_max = aoa_search_scope(2);
    aoa_display_step = aoa_search_scope(4);
    % ------ 2D ------
    figure;
    ax = pcolor(dist_search,aoa_search,aoa_dist_map);
    if display_flag_gt
        % ------ Ground Truth -----
        hold on;
        point = scatter(gt_mp_dists(1),gt_mp_aoas(1),80,'+r','linewidth',3);
        for mp_idx=2:num_of_mps
            hold on;
            scatter(gt_mp_dists(mp_idx),gt_mp_aoas(mp_idx),80,'+r','linewidth',3);
        end
        legend(point,'Groundtruth');
    end
    set(ax, 'LineStyle','none');
    xlabel('Range (m)');
    xlim([dist_min dist_max]);
    xticks(dist_min:dist_display_step:dist_max);
    ylabel('AoA (Deg)');
    ylim([aoa_min aoa_max]);
    yticks(aoa_min:aoa_display_step:aoa_max);
    title('Super-resolution Range-AoA');
    set(gca,'linewidth',1.5,'fontsize',20,'fontname','Arial');
    shading interp;
    colorbar;
    
    % ------ 3D ------
%     figure;
%     ax = surf(dist_search,aoa_search,aoa_dist_map);
%     if display_flag_gt
%         % ------ Ground Truth -----
%         hold on;
%         point = scatter3(gt_mp_dists(1),gt_mp_aoas(1),peak_amp_vec(peak_amp_idx(1)),'o','MarkerEdgeColor',[0.8500, 0.3250, 0.0980],'markerfacecolor','w','linewidth',2);
%         for mp_idx=2:num_of_mps
%             hold on;
%             scatter3(gt_mp_dists(mp_idx),gt_mp_aoas(mp_idx),peak_amp_vec(peak_amp_idx(mp_idx)),'o','MarkerEdgeColor',[0.8500, 0.3250, 0.0980],'markerfacecolor','w','linewidth',2);
%         end
%         legend(point,'Groundtruth','Location','best');
%     end
%     set(ax, 'LineStyle','none');
%     xlabel('Range (m)');
%     xlim([dist_min dist_max]);
%     xticks(dist_min:dist_display_step:dist_max);
%     ylabel('AoA (Deg)');
%     ylim([aoa_min aoa_max]);
%     yticks(aoa_min:aoa_display_step:aoa_max);
%     zlabel('Amplitude');
%     title('Super-resolution Range-AoA');
%     set(gca,'linewidth',1.5,'fontsize',20,'fontname','Arial');
%     colormap(jet);
%     shading interp;
%     colorbar;
end
end

