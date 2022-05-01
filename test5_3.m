N_fft=length(r)+length(mult_channel)-1-90;
r_fft=fft(r,N_fft);
mult_channel_fft=fft(mult_channel,N_fft);
r_circular=ifft(r_fft.*mult_channel_fft);

min_length=min(length(r_circular),length(r_linear));
err=r_circular(1:min_length)-r_linear(1:min_length);
err_am=abs(err);
figure;
plot(err_am(1:100));
max_err=max(abs(err));