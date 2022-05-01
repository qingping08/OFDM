%对称
clear;
clc;
carrier_count = 200; % 子载波数
symbol_count = 1;
ifft_length = 512;
CP_length = 128;
bit_per_symbol = 2; % qpsk调制

bit_length = carrier_count*symbol_count*bit_per_symbol;
sequence=randi([0 1],1,bit_length);
x_n= 2*sequence(1:2:end)+sequence(2:2:end);
const=[1+1i,-1+1i,1-1i,-1-1i];
qpsk=genqammod(x_n,const);
x_length=length(qpsk);

signal_time = ifft(qpsk);
% y_dft=fft([zeros(1,ifft_length-carrier_count),qpsk],ifft_length);%添加fft模块
y_dft=fft([qpsk],ifft_length);%添加fft模块
signal_time_dft=ifft(y_dft);
%==============加CP==================%
signal_time_CP = [signal_time(end-CP_length+1:end),signal_time];
signal_time_CP_dft = [signal_time_dft(end-CP_length+1:end),signal_time_dft];

y=reshape(signal_time_CP,1,[]);%传送信号
y_dft=reshape(signal_time_CP_dft,1,[]);%传送信号 dft

y_mag=abs(y.^2);
y_mag_dft=abs(y_dft.^2);
%===========计算峰均比=================%
% y_ifft_window=abs(y_window.^2);
figure;
subplot(2,1,1);
plot(y_mag);
mp=mean(y_mag);
pp=max(y_mag);
papr=10*log10(pp/mp);
title(['PAPR=',num2str(papr)]);
xlabel("不加DFT");
ylabel("信号幅度");
% ylim([0 max(y_mag)]);

% figure;
subplot(2,1,2);
plot(y_mag_dft);
mp1=mean(y_mag_dft);
pp1=max(y_mag_dft);
papr1=10*log10(pp1/mp1);
title(['PAPR=',num2str(papr1)]);
xlabel("加DFT");
ylabel("信号幅度");
ylim([0 2.1]);

%===========经过信道=============%
%信号保持
% t=0:1/Fs:1.2*T;
% signal=0;
% for i=1:N
%     signal_sub=rectpuls(t-(i-1/2)*T_a,T_a)*y(i);
%     signal=signal+signal_sub;
% end
% s_ifft=abs(signal.^2);
% figure;
% plot(t,s_ifft);
% mps=mean(s_ifft);
% pps=max(s_ifft);
% paprs=10*log10(pps/mps);
% title(num2str(paprs));
%====================================%
% function window=rcoswindow(alpha,bit_length)
%     warning off;
%     window = zeros(1,bit_length/2);
%     t = 1:bit_length/2;
%     T = bit_length/(2*(1+alpha));
%     window(t) = 0.5*(1 - sin(pi/(2*alpha*T)*(t-T)));
%     window(1:(1-alpha)*T) = 1;
%     window=[fliplr(window) window]';
% end