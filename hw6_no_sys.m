%非对称
clear;
clc;
carrier_count = 512; % 子载波数
symbol_count = 10;
ifft_length = 512;
CP_length = 128;
bit_per_symbol = 2; % qpsk调制

bit_length = carrier_count*symbol_count*bit_per_symbol;
sequence=randi([0 1],1,bit_length);
x_n= 2*sequence(1:2:end)+sequence(2:2:end);
const=[1+1i,-1+1i,1-1i,-1-1i];
qpsk=genqammod(x_n,const);
qpsk=qpsk.';%列向量 非共轭转置
x_length=length(qpsk);

bit_moded = reshape(qpsk,carrier_count,symbol_count);%串-并

signal_time = ifft(bit_moded,ifft_length);
%加CP
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

%===========计算峰均比=================%
y_ifft=abs(y.^2);
figure;
plot(y_ifft);