clear;
clc;
N=2048;
sequence=randi([0 1],1,N);
qpsk= 2*sequence(1:2:end)+sequence(2:2:end);
const=[1+1i,-1+1i,1-1i,-1-1i];
qpsk=genqammod(qpsk,const);
N=length(qpsk);%1024

B=3000;
delta_f=B/N;
T=1/delta_f;
Fs=30000;
T_a=1/B;%信号间隔
%x(n)
y=ifft(qpsk);
% plot(abs(y));
% xlim([1,1024])
% title('IFFT幅度')

%信号保持
t=0:1/Fs:T-1/Fs;
signal=0;
for i=1:N
    %signal_sub=rect((i-1)*T_a,i*T_a,t)*y(i);
    signal_sub=rectpuls(t-(i-1/2)*T_a,T_a)*y(i);
    signal=signal+signal_sub;
end
figure;
plot(t,abs(signal));
title('IFFT幅度')
figure;
h=myfilter2;
ylp=filter(h,signal);
plot(t,abs(ylp));
title('经过低通滤波器后')

figure;
plot(t,abs(signal),'r');
hold on;
plot(t,abs(ylp),'k');
legend('原信号','低通滤波')
xlim([0.22,0.23]);
title('取0.22-0.23段观察')

figure;
plot(t,abs(signal),'r');
hold on;
plot(t-T_a,abs(ylp),'k');
legend('原信号','低通滤波')
xlim([0.22,0.23]);
title('补偿T_a后')