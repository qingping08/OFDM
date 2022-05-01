clear;
clc;
N=20480;
sequence=randi([0 1],1,N);
qpsk= 2*sequence(1:2:end)+sequence(2:2:end);
x_n=qpsk;
const=[1+1i,-1+1i,1-1i,-1-1i];
qpsk=genqammod(qpsk,const);
N=length(qpsk);%1024

B=3000;
delta_f=B/N;
T=1/delta_f;
Fs=30000;
T_a=1/B;%信号间隔
f_c=15000;%调制频率
SNR=10;%信噪比
%x(n)
y=ifft(qpsk);
% plot(abs(y));
% xlim([1,1024])
% title('IFFT幅度')

%信号保持
t=0:1/Fs:1.2*T;
signal=0;
for i=1:N
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
xlim([0.22,0.23]);
title_str=['通过信噪比',num2str(SNR),'dB的信道'];
title(title_str);

r_low_t=filter(h,r_t);
figure;
plot(t,abs(s_t),'r');
hold on;
plot(t-T_a,abs(r_low_t),'k');
legend('通过信道前','通过接收端低通滤波器后')
xlim([0.22,0.23]);
title_str=['通过信噪比',num2str(SNR),'dB的信道'];
title(title_str);

figure;
plot(t,abs(signal),'g');
hold on;
plot(t,abs(s_t),'r');
hold on;
plot(t,abs(r_low_t),'k');
legend('IFFT序列','通过发送端低通滤波器后','通过接收端低通滤波器后')
xlim([0.22,0.23]);
title_str=['通过信噪比',num2str(SNR),'dB的信道'];
title(title_str);
%采样
r=1:N;
for i=1:N
 r(i)=r_low_t(fix((i-1/2+2)*T_a*Fs));
end
%FFT
%解调
Y=fft(r);
y_p=genqamdemod(Y,const);
figure;
scatter(real(Y),imag(Y),5,'filled');
title_str=['通过信噪比',num2str(SNR),'dB的信道'];
title(title_str);
fprintf(num2str(length(find(x_n==y_p))));