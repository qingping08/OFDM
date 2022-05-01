function [Y_sig,rate] = channel_f_diff(x_n,epsilon)
%对称
carrier_count = 200; % 子载波数
symbol_count = 10;
ifft_length = 512;
%bit_per_symbol = 2; % qpsk调制
CP_length = 128;

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

f_diff=epsilon*delta_f;
f_c_diff=f_c-f_diff;
r_t=r_t2.*exp(-1*1i*2*pi*f_c_diff*t);

r_low_t=filter(h,r_t);

%采样
r=1:N;
for i=1:N
 r(i)=r_low_t(fix((i-1/2+2)*T_a*Fs));
end

%==============串并转换=============%
Rx_data_sig = reshape(r,ifft_length+CP_length,[]);
%==========去CP==================%
Rx_data_sig(1:CP_length,:) = [];

%==============FFT解调===========%
Y_sig=fft(Rx_data_sig,ifft_length);

%========降采样===============%
data_sig = Y_sig(carrier_position,:);

Y_sig=reshape(data_sig,[],1).';%非共轭转置
%==========解调===============%
y_p_sig=genqamdemod(Y_sig,const);
ber_sig=1-length(find(x_n==y_p_sig))/x_length;

rate=ber_sig;
end

