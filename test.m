clear;
clc;
B=3000;
fft_n=4096;
sequence=randi([0 1],1,fft_n);
qpsk= 2*sequence(1:2:end)+sequence(2:2:end);
const=[1+1i,-1+1i,1-1i,-1-1i];
qpsk=genqammod(qpsk,const);
y=ifft(qpsk);
fft_n=length(qpsk);
delta_f=B/fft_n;
n=1:fft_n;
f=delta_f*(n-1);
T=1/delta_f;
fs=8000;
t=(-fft_n/2)/fs:1/fs:(fft_n/2-1)/fs;
signal_total=0;
figure;
for i=1:fft_n
    signal=rectpuls(t-(1/2)*T,T).*exp(1i*2*pi*f(i)*t);%*y(i);
    signal_total=signal_total+signal;
end
xdft = fft(signal_total);
xdft = xdft(1:fft_n/2+1);
psdx = (1/(fs*fft_n)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:fs/fft_n:fs/2;

plot(freq,10*log10(psdx));
hold on;

fft_n=4096*2;
sequence=randi([0 1],1,fft_n);
qpsk= 2*sequence(1:2:end)+sequence(2:2:end);
const=[1+1i,-1+1i,1-1i,-1-1i];
qpsk=genqammod(qpsk,const);
y=ifft(qpsk);
fft_n=length(qpsk);
delta_f=B/fft_n;
n=1:fft_n;
f=delta_f*(n-1);
T=1/delta_f;
fs=8000;
t=(-fft_n/2)/fs:1/fs:(fft_n/2-1)/fs;
signal_total=0;
for i=1:fft_n
    signal=rectpuls(t-(1/2)*T,T).*exp(1i*2*pi*f(i)*t);%*y(i);
    signal_total=signal_total+signal;
end
xdft = fft(signal_total);
xdft = xdft(1:fft_n/2+1);
psdx = (1/(fs*fft_n)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:fs/fft_n:fs/2;

plot(freq,10*log10(psdx));
hold on;
xlim([0,3500]);
legend("2048","4096");
title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
