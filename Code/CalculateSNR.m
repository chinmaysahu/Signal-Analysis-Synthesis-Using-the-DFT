% function that calculates SNR given original audio, estimated audio and framelength N=(64,128,256) 
function [Overall_SNR_1,Segmental_SNR_1,Overall_SNR_2,Segmental_SNR_2]=CalculateSNR(audio_signal,estimated_audio_signal_1,estimated_audio_signal_2,N)

Snr_I_Num=0;
Snr_I_Den=0;
Snr_Seg_I=0;

Snr_II_Num=0;
Snr_II_Den=0;
Snr_Seg_II=0;

L=length(audio_signal)/N;

for j=1:N:length(audio_signal) % 1:64:32768
  
  Xt=audio_signal(j:j+N-1); % Lth frame 
  estimatedXt_1=estimated_audio_signal_1(j:j+N-1);% estimated Xt method_1
  estimatedXt_2=estimated_audio_signal_2(j:j+N-1);% estimated Xt method_2
 
% Method 1
et=(Xt-estimatedXt_1);% Lth frame - estimated Lth frame
Snr_I_Num=Snr_I_Num+(Xt.'*Xt);
Snr_I_Den=Snr_I_Den+(et.'*et);

Snr_Seg_I=Snr_Seg_I+(10/L)*log10((Xt.'*Xt)/(et.'*et));  

%Method-2
et=(Xt-estimatedXt_2); 

Snr_II_Num=Snr_II_Num+(Xt.'*Xt);
Snr_II_Den=Snr_II_Den+(et.'*et);
Snr_Seg_II=Snr_Seg_II+((10/L)*log10(((Xt).'*Xt)/(et.'*et)));

end    

% Overall SNR
Overall_SNR_1=10*log10(Snr_I_Num/Snr_I_Den);
Overall_SNR_2=10*log10(Snr_II_Num/Snr_II_Den);

% Segmental SNR Dom Pick
Segmental_SNR_1=Snr_Seg_I;
% Segmental SNR Method-2
Segmental_SNR_2=Snr_Seg_II;

end