%对称
%频偏
clear;
clc;
bit_per_symbol = 2; % qpsk调制
carrier_count = 200; % 子载波数
symbol_count = 10;
ifft_length = 512;
bit_length = carrier_count*symbol_count*bit_per_symbol;
sequence=randi([0 1],1,bit_length);
x_n= 2*sequence(1:2:end)+sequence(2:2:end);

epsilon=1:30;
epsilon=epsilon*0.005;

ber_rate=zeros(1,length(epsilon));
for i=1:length(epsilon)
    [Y_sig,rate] = channel_f_diff(x_n,epsilon(i));
    ber_rate(i)=rate;
end
figure;
plot(epsilon,ber_rate);
xlabel("归一化频率");
ylabel("误码率");
title("误码率受频偏影响");