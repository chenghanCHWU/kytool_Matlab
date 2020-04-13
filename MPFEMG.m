function [MPF,EMG_s,alpha,beta,theta,delta]=MPFEMG(EEG,EMG,hour)
samplerate_EEG=125;
samplerate_EMG=250;

window=16;

minute=60;
second=60;


fo=2;
fc=32;

window_length_EEG=samplerate_EEG*window;
window_length_EMG=samplerate_EMG*window;
fft_E=[];fft_M=[];
slide_time=floor(hour*minute*second/window)*2-1;

for i=1:slide_time
     if i==1
        EEGsignal=double(EEG((i-1)*window_length_EEG+1:i*window_length_EEG));
        EMGsignal=double(EMG((i-1)*window_length_EMG+1:i*window_length_EMG)); 
     else
        EEGsignal=double(EEG((i-1)*window_length_EEG/2+1:(i+1)*window_length_EEG/2));
        EMGsignal=double(EMG((i-1)*window_length_EMG/2+1:(i+1)*window_length_EMG/2)); 
     end    
    

    [fM,EEG_fft]=fftpower(samplerate_EEG,EEGsignal);
    [fG,EMG_fft]=fftpower(samplerate_EMG,EMGsignal);
    
    fft_E=EEG_fft';
    fft_M=EMG_fft';
    
    
    if i==1
     
        f_index=[];
        for j=1 :round(samplerate_EEG/2)-1
            f_index_temp=find(fM>=(j-1) & fM<j);
            f_index=[f_index;f_index_temp];
        end
        
    end
    
    for j=1:round(samplerate_EEG/2)-1
        power_E(i,j)=sum(fft_E(f_index(j,:)));
        power_M(i,j)=sum(fft_M(f_index(j,:)));
    end
    
    EMG_s(i)=sum(power_M(i,32:62));
    MPF_temp=0;
    for j=fo:fc%round(samplerate_EEG/2)-1
        MPF_temp=MPF_temp+j*power_E(i,j);
    end
    MPF(i)=MPF_temp/sum(power_E(i,fo:fc));
    alpha(i)=sum(power_E(i,10:13));
    beta(i)=sum(power_E(i,13:32));
    theta(i)=sum(power_E(i,6:10));
    delta(i)=sum(power_E(i,1:4));
    
end

end

%beta 13-32 alpha 10-13 theta 6-10 MPF 2-32 delta 0.5-4 EMG 32-62