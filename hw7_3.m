clear;
clc;
f_c=300;
B=300;
N=2;
delta_f=B/N;
fft_n=4096*2;
n=1:N;
f=delta_f*(n-1);
f(2)=f(2)-delta_f*0.05;
T=1/delta_f;
fs=3000;
t=-fft_n/(2*fs):1/fs:(fft_n-1)/(2*fs);

%signal_total=0;
figure;
f_i=(0:fft_n-1)*fs/fft_n;
for i=1:N
    signal=rectpuls(t-(1/2)*T,T).*exp(1i*2*pi*f(i)*t);
    %signal_total=signal_total+signal;
    y=abs(fft(signal));
    %f_i=(0:fft_n-1)*fs/fft_n;
    plot(f_i,y)
    xlim([0,3.5*B]);
    %ylim([0,300]);
    hold on;
end
legend("f1","f2-f-diff");
% plot(f_i,abs(fft(signal_total)));
% xlim([0,B+50]);
% title("N=50")

