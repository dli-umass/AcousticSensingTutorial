function [rece_sw, rece_freq] = generate_receive_sw(amp,init_phase,vel,prop_delay,paras)
%GENERATEFMCWSOUNDWAVE Generate receive soundwave
%   amp             the attenuation factor
%   init_phase      initial phase 
%   vel             velocity
%   prop_delay      propagation delay in the beginning
%   paras           parameters 
%% Parameters
Fs = paras.fmcw_config.Fs;
Fc = paras.fmcw_config.Fc;
B = paras.fmcw_config.B;
T = paras.fmcw_config.T;  
Vs = paras.fmcw_config.Vs;
total_time = paras.fmcw_config.total_time;
%% Generate receive soundwave
num_of_samples = total_time*Fs;
rece_sw = zeros(1,num_of_samples);
rece_freq = zeros(1,num_of_samples);
for sam_idx=0:num_of_samples-1
    delay = prop_delay+2*vel*sam_idx/Fs/Vs;  
    if sam_idx/Fs >= delay
        t_n = mod(sam_idx/Fs-delay,T);    
        rece_freq(sam_idx+1) = Fc+B/T*t_n;
        phase = 2*pi*Fc*t_n+pi*B/T*t_n^2+init_phase;
        rece_sw(sam_idx+1) = amp*cos(phase);
    end
end
end