function [trans_sw_cos, trans_sw_sin, trans_freq] = generate_transmit_sw(A,init_phase,paras)
%GENERATEFMCWSOUNDWAVE Generate transmit sound or receive sound
%   A  amplitude
%   init_phase initial phase
%   Fs sampling rate
%   Fc start sweep frequency
%   B  bandwith
%   T  sweep time
%   total_time length of the sound wave 
%% Parameters
Fs = paras.fmcw_config.Fs;
Fc = paras.fmcw_config.Fc;
B = paras.fmcw_config.B;
T = paras.fmcw_config.T;  
%% Generate transmit soundwave
t = 0:1/Fs:T-1/Fs;
    
% Self-defined linear FMCW wave
trans_freq = Fc+B/T*t;
phase = Fc*t+B*power(t,2)/(2*T)+init_phase;
trans_sw_cos = A*cos(2*pi*phase).';
trans_sw_sin = A*sin(2*pi*phase).';
end