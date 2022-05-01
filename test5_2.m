% ==================================================
% 作者： 肆拾伍
% 时间：11/30 2019
% 版本：V3
% =================OFDM仿真参数说明：================
% 
% 子载波数 carrier_count     ---200
% 总符号数 symbol_count      ---100
% IFFT长度 ifft_length       ---512
% 循环前缀 CP_length         ---512/4=128
% 循环后缀 CS_length         ---20
% 升余弦窗系数 alpha         ---7/32
% 调制方式                   QAM16、QPSK 可选
% 多径幅度 mult_path_am      ---[1 0.2 0.1] 
% 多径时延 mutt_path_time    ---[0 20 50]
% ====================仿真过程=======================
% 产生0-1随机序列 => 串并转换 => 映射 => 取共轭、过采样
% => IFFT => 加循环前缀和后缀 => 加窗 => 并串转换 => 
% 多径信道 =>  加AWGN => 串并转换 => 去前缀 => FFT =>
% 下采样 =>  解映射 => 求误码率 
% ==================================================
clear all;
close all;
carrier_count = 200; % 子载波数
symbol_count = 100;
ifft_length = 512;
CP_length = 128;
CS_length = 20;
rate = [];
SNR =20;
bit_per_symbol = 4; % 调制方式决定
alpha = 1.5/32;
% ================产生随机序列=======================
bit_length = carrier_count*symbol_count*bit_per_symbol;
bit_sequence = round(rand(1,bit_length))'; % 列向量
% ================子载波调制方式1========================
% 1-28置零 29-228有效 229-285置零 286-485共轭 486-512置零
carrier_position = 29:228;
conj_position = 485:-1:286;
bit_moded = qammod(bit_sequence,16,'InputType','bit');
% ===================IFFT===========================
% =================串并转换==========================
ifft_position = zeros(ifft_length,symbol_count);
bit_moded = reshape(bit_moded,carrier_count,symbol_count);
ifft_position(carrier_position,:)=bit_moded(:,:);
ifft_position(conj_position,:)=conj(bit_moded(:,:));
signal_time = ifft(ifft_position,ifft_length);
% ==================加循环前缀和后缀==================
signal_time_C = signal_time;
% =======================加窗========================
signal_window = zeros(size(signal_time_C));
% 通过矩阵点乘
signal_window = signal_time_C.*repmat(rcoswindow(alpha,size(signal_time_C,1)),1,symbol_count);
% ===================发送信号，多径信道====================
signal_Tx = reshape(signal_window,1,[]); % 变成时域一个完整信号，待传输
signal_origin = reshape(signal_time_C,1,[]); % 未加窗完整信号
mult_path_am = [1 0.2 0.1]; %  多径幅度
mutt_path_time = [0 20 50]; % 多径时延
windowed_Tx = zeros(size(signal_Tx));
path2 = 0.2*[zeros(1,20) signal_Tx(1:end-20) ];
path3 = 0.1*[zeros(1,50) signal_Tx(1:end-50) ];
signal_Tx_mult = signal_Tx + path2 + path3; % 多径信号
% =====================发送信号频谱========================
% ====================未加窗信号频谱=======================
% 每个符号求频谱再平均，功率取对数
% ====================加窗信号频谱=========================
orgin_aver_power = 20*log10(mean(abs(fft(signal_window'))));
% ========================加AWGN==========================
signal_power_sig = var(signal_Tx); % 单径发送信号功率
signal_power_mut = var(signal_Tx_mult); % 多径发送信号功率
SNR_linear = 10^(SNR/10);
noise_power_mut = signal_power_mut/SNR_linear;
noise_power_sig = signal_power_sig/SNR_linear;
noise_sig = randn(size(signal_Tx))*sqrt(noise_power_sig);
noise_mut = randn(size(signal_Tx_mult))*sqrt(noise_power_mut);
% noise_sig=0;
% noise_mut=0;
Rx_data_sig = signal_Tx+noise_sig;
Rx_data_mut = signal_Tx_mult+noise_mut;
% =======================串并转换==========================
Rx_data_mut = reshape(Rx_data_mut,ifft_length,[]);
Rx_data_sig = reshape(Rx_data_sig,ifft_length,[]);
% ====================去循环前缀和后缀======================
% =========================FFT=============================
fft_sig = fft(Rx_data_sig);
fft_mut = fft(Rx_data_mut);
% =========================降采样===========================
data_sig = fft_sig(carrier_position,:);
data_mut = fft_mut(carrier_position,:);
figure
scatter(real(reshape(data_sig,1,[])),imag(reshape(data_sig,1,[])),'.')
grid on;
figure
scatter(real(reshape(data_mut,1,[])),imag(reshape(data_mut,1,[])),'.')
grid on;
% =========================逆映射===========================
bit_demod_sig = reshape(qamdemod(data_sig,16,'OutputType','bit'),[],1);
bit_demod_mut = reshape(qamdemod(data_mut,16,'OutputType','bit'),[],1);
% =========================误码率===========================
error_bit_sig = sum(bit_demod_sig~=bit_sequence);
error_bit_mut = sum(bit_demod_mut~=bit_sequence);
error_rate_sig = error_bit_sig/bit_length;
error_rate_mut = error_bit_mut/bit_length;
rate = [rate; error_rate_sig error_rate_mut];
% ==========================================================
% ==========================================================
function window=rcoswindow(alpha,bit_length)
    warning off;
    window = zeros(1,bit_length/2);
    t = 1:bit_length/2;
    T = bit_length/(2*(1+alpha));
    window(t) = 0.5*(1 - sin(pi/(2*alpha*T)*(t-T)));
    window(1:(1-alpha)*T) = 1;
    window=[fliplr(window) window]';
end
