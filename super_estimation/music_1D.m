function [fom_freq,fom_sp] = music_1D(sig,Fs,freq_min,freq_max,freq_step,K,W)
%MUSIC_1D Estimate the frequency for one dimension signal
%   sig         input signals
%   Fs          sampling rate
%   freq_min    minimum frequency
%   freq_max    maximum frequency
%   freq_step   frequency stepsize
%   K           subsampling factor
%   W           length of smoothing window
%
%   This algorithm combines the subsampling and data smoothing:
%   1. Set K=1 to avoid subsampling 
%   2. Set W as the number of samples in the subarray after subsampling to avoid data smoothing 
%   Reference:
%   Subsampling 
%   [1] Indoor Follow Me Drone
%   Data Smoothing
%   [2] Superresolution Techniques for Time-Domain Measurements with a Network Analyzer
%   [3] Super-Resolution TOA Estimation With Diversity for Indoor Geolocation
%% Construct the correlation matrix
freq_num = ceil((freq_max-freq_min)/freq_step);           % number of estimated frequencies
N  = length(sig);                                         % number of samples
SUB_N = floor((N-K)/K)+1;                                 % number of samples after subsampling
NUM_W = SUB_N - W + 1;                                    % number of the smoothing window
trans_matrix = eye(W,W);
J = fliplr(trans_matrix);

X = zeros(W, NUM_W);
if SUB_N < NUM_W
    fprintf('Warning: The subsampling factor K is too large!\n');
    return;
end

R_sub = 0;                                                 % autocorrelation matrix after subsampling
for sub_idx=1:K
    % ------ subsampling ------
    sub_sig = sig(sub_idx:K:sub_idx+K*(SUB_N-1));
    
    % ------ data smoothing -----
    for win_idx=1:NUM_W
        X(:,win_idx) = sub_sig(win_idx:win_idx+W-1);
    end
    R_sub = R_sub + 1/(2*NUM_W)*(X*X'+J*(X*X').'*J);
end
R_sub = R_sub/K;

% ------ MUSIC -----
[V,D] = eig(R_sub);                     % Compute the eigenvalues and eigenvectors of the auto-correlation matrix
[D_val,I] = sort(diag(D),'descend');    % Sort the eigenvectors in a descending order in terms of the magnitude of corresponding eigenvalues.

% ----- Determine the noise number ------
Level=10^4;
if abs(max(D_val)/min(D_val))<Level
    fprintf('Warning: The eigen values may not correct!\n         The estimated frequicies may be wrong!\n')
end
noise_num = W - length(find(D_val<max(D_val)/Level));

% MDL_vec = zeros(W,1);
% for k=1:W
%     numerator = 1;
%     for i=k:W
%         numerator = numerator * D_val(i)^(1/(W+1-k));
%     end
%     denominator = 0;
%     for i=k:W
%         denominator = denominator + D_val(i);
%     end
%     denominator = (1/(W+1-k))*denominator;
%     MDL_vec(k) = -(NUM_W*(W+1-k))*log((numerator/denominator)) + 1/4*k*(2*(W+1)-k+1)*log(NUM_W);
% end
% [~,mp_num] = min(MDL_vec);
% noise_num = W - mp_num;

fom_freq = zeros(freq_num,1);  % Initialize the frequcy in frequency spectrum
fom_sp = zeros(freq_num,1);    % Initialize the frequency spectrum 
Ts = 1/Fs;
time = 0:K*Ts:(W-1)*K*Ts;      % Time vector for each sample in steering vector
for freq_idx=1:freq_num
    fom_freq(freq_idx) = freq_min + (freq_idx-1)*freq_step;
    s_f = exp(1i*2*pi*time*fom_freq(freq_idx)).';   % Steering vector
    Rn = V(:, I(noise_num+1:W));                   % Noise space matrix
    fom_sp(freq_idx) = 1/(s_f'*(Rn*Rn')*s_f);
end


end

