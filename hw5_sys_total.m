%对称
%同时包含加不加cp两种
clear;
clc;
carrier_count = 200; % 子载波数
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
signal_time_noCP = signal_time; %不加cp

y=reshape(signal_time_CP,1,[]);%传送信号
y_noCP=reshape(signal_time_noCP,1,[]);%传送信号 不加cp

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

%==========多径信道、单径信道==========%
mult_path_am = [1 0.2 0.3]; %  多径幅度
mult_path_time = [0 20 50]; % 多径时延
r_sig=r;%单径信道
r_mult=r;%多径信道
for i=2:length(mult_path_am)
    r_mult=r_mult+mult_path_am(i)*[zeros(1,mult_path_time(i)) r(1:end-mult_path_time(i)) ];
end

%==============串并转换=============%
Rx_data_mult = reshape(r_mult,ifft_length+CP_length,[]);
Rx_data_sig = reshape(r_sig,ifft_length+CP_length,[]);
Rx_data_sig_copy=Rx_data_sig;
%==========去CP==================%
Rx_data_sig(1:CP_length,:) = [];
Rx_data_mult(1:CP_length,:) = [];

%==============FFT解调===========%
Y_sig=fft(Rx_data_sig,ifft_length);
Y_mult=fft(Rx_data_mult,ifft_length);

%========降采样===============%
data_sig = Y_sig(carrier_position,:);
data_mult = Y_mult(carrier_position,:);

Y_sig_copy=Y_sig;
Y_sig=reshape(data_sig,[],1).';%非共轭转置
Y_mult=reshape(data_mult,[],1).';
%===========星座图===============%
figure;
subplot(2,1,1);
scatter(real(Y_sig),imag(Y_sig),5,'filled');
grid;
y_p_sig=genqamdemod(Y_sig,const);
ber_sig=1-length(find(x_n==y_p_sig))/x_length;
title_str=['单径信道 误码率为',num2str(ber_sig*100),'%'];
title(title_str);
xlim([-2,2]);
ylim([-2,2]);
subplot(2,1,2);
scatter(real(Y_mult),imag(Y_mult),5,'filled');
grid;
y_p_mult=genqamdemod(Y_mult,const);
ber_mult=1-length(find(x_n==y_p_mult))/x_length;
title_str=['多径信道 误码率为',num2str(ber_mult*100),'%'];
title(title_str);
xlim([-2,2]);
ylim([-2,2]);