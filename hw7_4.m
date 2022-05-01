clear;
clc;
symbol_count=1024;
bit_per_symbol=2;
bit_length=symbol_count*bit_per_symbol;
sequence=randi([0 1],1,bit_length);
x_n= 2*sequence(1:2:end)+sequence(2:2:end);
const=[1+1i,-1+1i,1-1i,-1-1i];
qpsk=genqammod(x_n,const);
N_fft=length(qpsk);
CFO=0:1024;
CFO=CFO*0.1;
ICI=zeros(1,length(CFO));
for i=1:length(CFO)
    x_ifft=ifft(qpsk,N_fft); %频域变时域
    y_CFO=add_CFO(x_ifft,CFO(i),N_fft);%加入频偏
    y_fft=fft(y_CFO);
    err=y_fft-qpsk;
    ICI(i)=abs(err(1))^2;
end
figure;
plot(CFO/1024,ICI);
xlabel("归一化频偏");
ylabel("ICI");
title("ICI随频偏变化")
%实现对OFDM加入频偏
function y_CFO=add_CFO(y,CFO,Nfft)
% nn=0:length(y)-1;
% y_CFO=y.*exp(1i*2*pi*CFO*nn/Nfft);
y_CFO=y;
y_CFO(1)=y_CFO(1)*exp(-1i*2*pi*1*CFO/Nfft);
end