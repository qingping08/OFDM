clear;
clc;
N=512*10*2;
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
y=ifft(qpsk);

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
mult_path_am = [1 0.2 0.4]; %  多径幅度
mult_path_time = [0 20 50]; % 多径时延
r_mult=r;
for i=2:length(mult_path_am)
    r_mult=r_mult+mult_path_am(i)*[zeros(1,mult_path_time(i)) r(1:end-mult_path_time(i)) ];
end
%FFT
%解调
Y=fft(r);
Y_mult=fft(r_mult);
figure;
scatter(real(Y),imag(Y),5,'filled');
grid;
y_p=genqamdemod(Y,const);
ber_sig=1-length(find(x_n==y_p))/N;
title_str=['单径信道 误码率为',num2str(ber_sig*100),'%'];
title(title_str);
xlim([-2,2]);
ylim([-2,2]);
figure;
scatter(real(Y_mult),imag(Y_mult),5,'filled');
grid;
y_p_mult=genqamdemod(Y_mult,const);
ber_mult=1-length(find(x_n==y_p_mult))/N;
title_str=['多径信道 误码率为',num2str(ber_mult*100),'%'];
title(title_str);
xlim([-2,2]);
ylim([-2,2]);
%fprintf(num2str(length(find(x_n==y_p_mult))));