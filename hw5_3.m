clc;
clear;
carrier_count = 200; % 子载波数
symbol_count = 10;
ifft_length = 512;
bit_per_symbol = 2; % qpsk调制

bit_length = carrier_count*symbol_count*bit_per_symbol;
sequence=randi([0 1],1,bit_length);
x_n= 2*sequence(1:2:end)+sequence(2:2:end);

N=10;
CP_step=10;
CP_length=1:N;
CP_length=CP_length*CP_step;

ber_mult=zeros(1,N);
for i=1:N
    [Y_sig,Y_mult,rate] = channel(CP_length(i),x_n);
    ber_mult(i)=rate(2);
end
figure;
plot(CP_length,ber_mult*100);
xlabel("CP length");
ylabel("ber rate/%");
title("误码率与CP长度的关系");