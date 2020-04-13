function [f_length,power_amp]=fftpower(samplerate,sig)

% vaild_frequency=round(samplerate/2);% window size 64 second
datalength = length(sig);
index=1:datalength/2;
hammingwindow = hamming(datalength);                                    % hamming window
f_length = samplerate*linspace(0,1,datalength);                                 % frequency calculate   
s_temp=detrend(sig);
s_temp=s_temp .* hammingwindow;
sfft = fft(s_temp,datalength);
power_amp=abs(sfft(index));
theta=mod(angle(sfft(index))*180/pi,360);
f_length=f_length(index);

end