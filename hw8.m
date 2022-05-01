clear;
clc;
symbol_count=200;
bit_per_symbol=2;
N=symbol_count*bit_per_symbol;
sequence=randi([0 1],1,N);
x_n= 2*sequence(1:2:end)+sequence(2:2:end);
const=[1+1i,-1+1i,1-1i,-1-1i];
qpsk=genqammod(x_n,const);
N=length(qpsk);%1024

ifft_length=512;
B=3000;
delta_f=B/ifft_length;
T=1/delta_f;
Fs=30000;
T_a=1/B;%信号间隔
f_c=15000;%调制频率
SNR=10;%信噪比
%x(n)
y_dft=fft(qpsk,ifft_length);%添加fft模块
y=ifft(y_dft);
% plot(abs(y));
% xlim([1,1024])
% title('IFFT幅度')

%信号保持
t=0:1/Fs:1.2*T;
signal=0;
for i=1:ifft_length
    %signal_sub=rect((i-1)*T_a,i*T_a,t)*y(i);
    signal_sub=rectpuls(t-(i-1/2)*T_a,T_a)*y(i);
    signal=signal+signal_sub;
end

h=myfilter2;
s_t=filter(h,signal);

s_t2=s_t.*exp(1i*2*pi*f_c*t);
r_t2=awgn(s_t2,SNR,'measured');

r_t=r_t2.*exp(-1*1i*2*pi*f_c*t);
figure;
plot(t-T_a,abs(s_t),'r');
hold on;
plot(t-T_a,abs(r_t),'k');
legend('通过信道前','通过信道后')
xlim([0.05,0.06]);
title_str=['通过信噪比',num2str(SNR),'dB的信道'];
title(title_str);

r_low_t=filter(h,r_t);
figure;
plot(t,abs(s_t),'r');
hold on;
plot(t-T_a,abs(r_low_t),'k');
legend('通过信道前','通过接收端低通滤波器后')
xlim([0.05,0.06]);
title_str=['通过信噪比',num2str(SNR),'dB的信道'];
title(title_str);

figure;
plot(t,abs(signal),'g');
hold on;
plot(t,abs(s_t),'r');
hold on;
plot(t,abs(r_low_t),'k');
legend('IFFT序列','通过发送端低通滤波器后','通过接收端低通滤波器后')
xlim([0.05,0.06]);
title_str=['通过信噪比',num2str(SNR),'dB的信道'];
title(title_str);
%采样
r=1:ifft_length;
for i=1:ifft_length
 r(i)=r_low_t(fix((i-1/2+2)*T_a*Fs));
end
%FFT
%解调
Y=fft(r);
Y_idft=ifft(Y,ifft_length);
Y_final=Y_idft(1:symbol_count);
y_p=genqamdemod(Y_final,const);
figure;
scatter(real(Y_final),imag(Y_final),5,'filled');
title_str=['通过信噪比',num2str(SNR),'dB的信道'];
title(title_str);
fprintf(num2str(length(find(x_n==y_p))));