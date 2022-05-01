%对称
clear;
clc;
carrier_count = 200; % 子载波数
symbol_count = 1;
ifft_length = 512;
CP_length = 128;
bit_per_symbol = 2; % qpsk调制

mult_path_am = [1 0.2 0.4]; %  多径幅度
mult_path_time = [0 20 50]; % 多径时延

bit_length = carrier_count*symbol_count*bit_per_symbol;
sequence=randi([0 1],1,bit_length);
x_n= 2*sequence(1:2:end)+sequence(2:2:end);
const=[1+1i,-1+1i,1-1i,-1-1i];
qpsk=genqammod(x_n,const);
qpsk=qpsk.';%列向量 非共轭转置
x_length=length(qpsk);

%==========串并转换===============%
bit_moded = reshape(qpsk,carrier_count,symbol_count);%串-并
% 1-28置零 29-228有效 229-285置零 286-485共轭 486-512置零
carrier_position = 29:228;
conj_position = 485:-1:286;
ifft_position = zeros(ifft_length,symbol_count);
ifft_position(carrier_position,:)=bit_moded(:,:);
ifft_position(conj_position,:)=conj(bit_moded(:,:));
signal_time = ifft(ifft_position,ifft_length);

%==============加CP==================%
signal_time_CP = [signal_time(end-CP_length+1:end,:);signal_time];

y=reshape(signal_time_CP,1,[]);%传送信号
y_mag=abs(y.^2);
%=============传送参数===============%
N=length(y);
B=3000;
delta_f=B/N;
T=1/delta_f;
Fs=30000;%采样率
T_a=1/B;%信号间隔
f_c=15000;%调制频率  
SNR=10;%信噪比

%============加窗======================%
%signal_window = zeros(size(signal_time_CP));
% 通过矩阵点乘
% alpha=0.5;
% signal_window = signal_time_CP.*repmat(rcoswindow(alpha,size(signal_time_CP,1)),1,symbol_count);
% y_window=reshape(signal_window,1,[]);

%==============clipConv============%
y_max=0.7*max(y_mag);
y_mag1=zeros(1,N);
for j=1:N                                   %Clipping the signals above threshold(here 70% of original value)
if(y_mag(j)>y_max)
    y_mag1(j)=y_max;
else
    y_mag1(j)=y_mag(j);
end    
end

h=ones(1,N);
y_mag2=conv(y_mag1,h);
%===========计算峰均比=================%
% y_ifft_window=abs(y_window.^2);
figure;
subplot(3,1,1);
plot(y_mag);
mp=mean(y_mag);
pp=max(y_mag);
papr=10*log10(pp/mp);
title(['PAPR=',num2str(papr)]);
xlabel("传输信号(加CP)");
ylabel("信号幅度");
ylim([0 max(y_mag)]);

% figure;
subplot(3,1,2);
plot(y_mag1);
mp1=mean(y_mag1);
pp1=max(y_mag1);
papr1=10*log10(pp1/mp1);
title(['PAPR=',num2str(papr1)]);
xlabel("消波后");
ylabel("信号幅度");
ylim([0 max(y_mag)]);

% figure;
subplot(3,1,3);
plot(y_mag2);
mp2=mean(y_mag2);
pp2=max(y_mag2);
papr2=10*log10(pp2/mp2);
title(['PAPR=',num2str(papr2)]);
xlabel("卷积后");
ylabel("信号幅度");
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