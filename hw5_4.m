%线性卷积 循环卷积clear;
clc;
carrier_count = 200; % 子载波数
symbol_count = 10;
ifft_length = 512;
CP_length = 128;
bit_per_symbol = 2; % qpsk调制

% mult_path_am = [1 0.2 0.4]; %  多径幅度
% mult_path_time = [0 20 50]; % 多径时延

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

%传送参数
N=length(y);
B=3000;
delta_f=B/N;
T=1/delta_f;
Fs=30000;%采样率
T_a=1/B;%信号间隔
f_c=15000;%调制频率  
SNR=10;%信噪比

%===========经过信道=============%
%信号保持
t=0:1/Fs:1.2*T;
signal=0;
for i=1:N
    signal_sub=rectpuls(t-(i-1/2)*T_a,T_a)*y(i);
    signal=signal+signal_sub;
end

h=myfilter2;
s_t=filter(h,signal);

s_t2=s_t.*exp(1i*2*pi*f_c*t);
r_t2=awgn(s_t2,SNR,'measured');

r_t=r_t2.*exp(-1*1i*2*pi*f_c*t);

r_low_t=filter(h,r_t);

%采样
r=1:N;
for i=1:N
 r(i)=r_low_t(fix((i-1/2+2)*T_a*Fs));
end

%多径信道
mult_path_am = [1 0.2 0.4]; %  多径幅度
mult_path_time = [0 20 50]; % 多径时延

length_mult_channel=mult_path_time(end)+1;
mult_channel=zeros(1,length_mult_channel);
for i=1:length(mult_path_am)
    mult_channel(mult_path_time(i)+1)=mult_path_am(i);
end
%==========多径信道===============%
r_mult=r;
for i=2:length(mult_path_am)
    r_mult=r_mult+mult_path_am(i)*[zeros(1,mult_path_time(i)) r(1:end-mult_path_time(i)) ];
end
%===============线性卷积===================%
r_linear=conv(r,mult_channel);
%=============循环卷积=====================%
N_fft=length(r)+length(mult_channel)-1-60;
r_fft=fft(r,N_fft);
mult_channel_fft=fft(mult_channel,N_fft);
r_circular=ifft(r_fft.*mult_channel_fft);

min_length=min(length(r_circular),length(r_linear));
err=r_circular(1:min_length)-r_linear(1:min_length);
err_am=abs(err);
figure;
plot(err_am(1:100));
% plot(err_am);
ylabel("error");
title("线性卷积与循环卷积的误差")
max_err=max(abs(err));