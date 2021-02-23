function all_sig_path_paras = trad_aoa_FFT(sig,paras)
%TRAD_AOA_FFT Estimate the aoa using FFT
%   sig             input signal
%   paras           parameters

%% Parameters
Vs = paras.fmcw.Vs;
T = paras.fmcw.T;
Fc = paras.fmcw.Fc;
single_chirp_len = paras.fmcw.single_chirp_len;

num_of_chirps = paras.algos.num_of_chirps;
aoa_max = paras.algos.aoa_max;
aoa_fft_size = paras.algos.trad.aoa_fft_size;
display_flag = paras.algos.display_flag;

num_of_mps = paras.multipath.num_of_mps;
gt_mp_aoas = paras.multipath.aoas;
%% Algos
aoa_search = acosd(linspace(-1,1,aoa_fft_size));
aoa_idx = (aoa_search>=180-aoa_max) & (aoa_search<=aoa_max);
aoa_search = aoa_search(aoa_idx);

%% ------ aoa + samples -----
chirp_idx = 1;
aoa_samples = squeeze(sig(:,:,chirp_idx));
aoa_samples_fft = fftshift(fft(aoa_samples,aoa_fft_size,2),2);
aoa_samples_map = abs(aoa_samples_fft).^2/aoa_fft_size;

%% ------ aoa + chirps -----
sample_idx = 2000;
aoa_chirps = squeeze(sig(sample_idx,:,:));
aoa_chirps_fft = fftshift(fft(aoa_chirps,aoa_fft_size,1),1);
aoa_chirps_map = abs(aoa_chirps_fft).^2/aoa_fft_size;
%% Display
if display_flag
%     celldisp(all_sig_path_paras);
    
    % ------ 1D ------
    figure;
    plot(aoa_search,aoa_chirps_map(:,1),'linewidth',3);
    for mp_idx=1:num_of_mps
        hold on;
        plot([gt_mp_aoas(mp_idx) gt_mp_aoas(mp_idx)],[min(aoa_chirps_map(:,1)) max(aoa_chirps_map(:,1))],'--r','linewidth',3);
    end
    xlabel('Velocity (m/s)');
    ylabel('Amplitude');
    title('Velocity FFT');
    set(gca,'linewidth',1.5,'fontsize',20,'fontname','Times');
    
    % ------ 2D ------
    figure;
    ax = pcolor(1:single_chirp_len,aoa_search,aoa_samples_map.');
    for mp_idx=1:num_of_mps
        hold on;
        plot([1 single_chirp_len],[gt_mp_aoas(mp_idx) gt_mp_aoas(mp_idx)],'--w','linewidth',3);
    end
    set(ax, 'LineStyle','none');
    xlabel('Sample Number');
    ylabel('AoA (deg)');
    title('AoA-Samples FFT');
    set(gca,'linewidth',1.5,'fontsize',20,'fontname','Times');
    colormap(jet);
    shading interp;
    colorbar;
    
    % ------ 2D ------
    figure;
    ax = pcolor(1:num_of_chirps,aoa_search,aoa_chirps_map);
    for mp_idx=1:num_of_mps
        hold on;
        plot([1 num_of_chirps],[gt_mp_aoas(mp_idx) gt_mp_aoas(mp_idx)],'--w','linewidth',3);
    end
    set(ax, 'LineStyle','none');
    xlabel('Chirp Number');
    ylabel('AoA (deg)');
    title('AoA-Chirps FFT');
    set(gca,'linewidth',1.5,'fontsize',20,'fontname','Times');
    colormap(jet);
    shading interp;
    colorbar;
end
end

