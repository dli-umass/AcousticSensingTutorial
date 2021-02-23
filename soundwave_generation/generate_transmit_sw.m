function [trans_sw_cos, trans_sw_sin, trans_freq] = generate_transmit_sw(A,init_phase,Fs,Fc,B,T,total_time)
%GENERATEFMCWSOUNDWAVE Generate transmit sound or receive sound
%   A  amplitude
%   init_phase initial phase
%   Fs sampling rate
%   Fc start sweep frequency
%   B  bandwith
%   T  sweep time
%   total_time length of the sound wave 

t = 0:1/Fs:total_time-1/Fs;
t_n = mod(t, T);  % time for the nth chirp
    
% Self-defined linear FMCW wave
trans_freq = Fc+B/T*t_n;
phase = Fc*t_n+B*power(t_n,2)/(2*T)+init_phase;
trans_sw_cos = A*cos(2*pi*phase);
trans_sw_sin = A*sin(2*pi*phase);
end