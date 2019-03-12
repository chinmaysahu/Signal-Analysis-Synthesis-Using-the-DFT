% Author: Chinmay Sahu
% Class: EE 501 - Digital Signal Processing
% Project: Signal Analysis/Synthesis using the DFT
% email: sahuc@clarkson.edu
% Clarkson Univeristy
% Nov 2017; Last revision: 19-Nov-2017
clc
clear
% close all
filename = 'cleanspeech.wav';% Loading the sound file

Reconstructed_array=[]; % initializing reconstructed audio array for method-2
Dom_recon_array=[];  % initializing reconstructed audio array for method-1
first_n_picks=[]; %initializing array to keep first n picks for method-2
peak_keepers=[]; %initializing array to keep first n peaks for method-1


SNR_64_I=[];
SNR_64_II=[];
Snr_I=0;
Snr_II=0;
Snr_Seg_I=0;
Snr_Seg_II=0;

for fftPoints= 1 : 3
    
switch fftPoints
case 1
N=64;% Framerate
case 2
N=128;
case 3
N=256;
end


[y,Fs] = audioread('cleanspeech.wav');% read the sound file

%N=64; % decideing the  N=64,128,256
L=length(y)/N;
y1=reshape(y,[N,L]);% reshaping the matrix to N-point frames of Lth columns 
first_n_picks=zeros(1+(N/2),1);
 % 1+(N/2) independent components, only n independent components are retained for data reconstruction (n<1+(N/2))
for n=1:1:(1+(N/2)) % n varies from 1: 1+(N/2) example: for N=64, n=1:33
    
Reconstructed_array=[];
Dom_recon_array=[];
Snr_I_Num=0;
Snr_I_Den=0;
Snr_II_Num=0;
Snr_II_Den=0;

Snr_II=0;
Snr_Seg_I=0;
Snr_Seg_II=0;

for l=1:1:L % Lth Frame 
    
N_frame_FFT=fft(y1(:,l),N); % N frame FFT

% Method 1 dominant signal pick
[sortedValues,sortIndex] = sort(abs(N_frame_FFT(1:1+(N/2))),'descend');  %# Sort the values in descending order                                      
dominant_peaks = sortIndex(1:n);  %Careful of symmetry Z= noumber of dominant peaks
peak_keepers=zeros(1+(N/2),1);
for q=1:1:n
peak_keepers(dominant_peaks(q))=N_frame_FFT(dominant_peaks(q)); % Keeping dominant peaks at their indices
end


Dom_FlipSig=flipud(conj(peak_keepers)); % flip
Dom_Sym_Sig=[peak_keepers;Dom_FlipSig(2:N/2)]; % Joining flipped signal
Dom_N_frame_IFFT=ifft(Dom_Sym_Sig,N); % Takinf IFFT of symmetric signal
Dom_recon_array=[Dom_recon_array;Dom_N_frame_IFFT]; % reconstructing audio data 



Xt=y1(:,l); % Lth frame
et=(Xt-Dom_N_frame_IFFT);% Lth frame - estimated Lth frame

Snr_I_Num=Snr_I_Num+(Xt.'*Xt);
Snr_I_Den=Snr_I_Den+(et.'*et);

Snr_Seg_I=Snr_Seg_I+(10/L)*log10((Xt.'*Xt)/(et.'*et)); 

%Method-2 the first n components are retained 

first_n_picks=N_frame_FFT(1:n); % picking first n picks
first_n_picks=[first_n_picks;zeros((1+(N/2)-n),1)]; % padding zeros 
FlipSig=flipud(conj(first_n_picks)); % fliping after taking conjugate
Sym_Sig=[first_n_picks;FlipSig(2:N/2)]; % Joining to get symmtry
N_frame_IFFT=ifft(Sym_Sig,N); % IFFT
Reconstructed_array=[Reconstructed_array;N_frame_IFFT]; % reconstruction

Xt=y1(:,l);
et=(Xt-N_frame_IFFT); 

Snr_II_Num=Snr_II_Num+(Xt.'*Xt);
Snr_II_Den=Snr_II_Den+(et.'*et);
Snr_Seg_II=Snr_Seg_II+((10/L)*log10(((Xt).'*Xt)/(et.'*et)));

end


Reconstructed_audio_array_I{n}=Dom_recon_array; % recontructed audio for n method 1
Reconstructed_audio_array_II{n}=Reconstructed_array; % recontructed audio for n method 2


% Overall SNR
SNR_I(n)=10*log10(Snr_I_Num/Snr_I_Den);
SNR_II(n)=10*log10(Snr_II_Num/Snr_II_Den);

% Segmental SNR Dom Pick
SNR_SEGMENTAL_I(n)=Snr_Seg_I;
% Segmental SNR
SNR_SEGMENTAL_II(n)=Snr_Seg_II;
end % end of for loop where we move n using Method-1 or Method-2 to find SNR ,audios



switch fftPoints
case 1 % retrive audios and SNR for 64 frame
    
%Retrived Audio Signals
Compressed_Audio_Signals_64_1=Reconstructed_audio_array_I;
Compressed_Audio_Signals_64_2=Reconstructed_audio_array_II;   
% Overall SNR    
SNR_64_1=SNR_I;
SNR_64_2=SNR_II;
% Segmental
SNR_SEGMENTAL_64_1=SNR_SEGMENTAL_I;
SNR_SEGMENTAL_64_2=SNR_SEGMENTAL_II;

case 2
    
%Retrived Audio Signals
Compressed_Audio_Signals_128_1=Reconstructed_audio_array_I;
Compressed_Audio_Signals_128_2=Reconstructed_audio_array_II;  
% Overall SNR    
SNR_128_1=SNR_I;
SNR_128_2=SNR_II;
% Segmental
SNR_SEGMENTAL_128_1=SNR_SEGMENTAL_I;
SNR_SEGMENTAL_128_2=SNR_SEGMENTAL_II;

case 3
    
%Retrived Audio Signals
Compressed_Audio_Signals_256_1=Reconstructed_audio_array_I;
Compressed_Audio_Signals_256_2=Reconstructed_audio_array_II;      
% Overall SNR    
SNR_256_1=SNR_I;
SNR_256_2=SNR_II;
% Segmental
SNR_SEGMENTAL_256_1=SNR_SEGMENTAL_I;
SNR_SEGMENTAL_256_2=SNR_SEGMENTAL_II;
end

end



% Overall SNR method 1
figure
grid on
box on
grid minor
hold on
plot(100*((1:33)/64),SNR_64_1,':k','LineWidth',2);
hold on
plot(100*((1:65)/128),SNR_128_1,'--b','LineWidth',2);
hold on
plot(100*((1:129)/256),SNR_256_1,':r','LineWidth',2);
hold on
title('Method-1,Overall SNR vs % of components selected (n/N)','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
xlabel('% of components selected (n/N)','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
ylabel('Overall SNR(dB)','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','bold');

legend({'N=64','N=128','N=256'},'FontName','Times New Roman','FontSize',12,'FontWeight','bold');


% Segmental SNR method 1
figure
grid on
box on
grid minor
hold on
plot(100*((1:33)/64),SNR_SEGMENTAL_64_1,':k','LineWidth',2);
hold on
plot(100*((1:65)/128),SNR_SEGMENTAL_128_1,'--b','LineWidth',2);
hold on
plot(100*((1:129)/256),SNR_SEGMENTAL_256_1,':r','LineWidth',2);


title('Method-1,Segmental SNR vs % of components selected (n/N)','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
xlabel('% of components selected (n/N)','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
ylabel('Segmental SNR(dB)','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','bold');
legend({'N=64','N=128','N=256'},'FontName','Times New Roman','FontSize',12,'FontWeight','bold');

% Overall SNR method 2
figure
grid on
box on
grid minor
hold on
plot(100*((1:33)/64),SNR_64_2,':k','LineWidth',2);
hold on
plot(100*((1:65)/128),SNR_128_2,'--b','LineWidth',2);
hold on
plot(100*((1:129)/256),SNR_256_2,':r','LineWidth',2);


title('Method-2,Overall SNR vs % of components selected (n/N) ','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
xlabel('% of components selected (n/N)','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
ylabel('Overall SNR(dB)','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','bold');
legend({'N=64','N=128','N=256'},'FontName','Times New Roman','FontSize',12,'FontWeight','bold');


% Segmental SNR method 2
figure
grid on
box on
grid minor
hold on
plot(100*((1:33)/64),SNR_SEGMENTAL_64_2,':k','LineWidth',2);
hold on
plot(100*((1:65)/128),SNR_SEGMENTAL_128_2,'--b','LineWidth',2);
hold on
plot(100*((1:129)/256),SNR_SEGMENTAL_256_2,':r','LineWidth',2);

title('Method-2,Segmental SNR vs % of components selected (n/N)','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
xlabel('% of components selected (n/N)','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
ylabel('Segmental SNR(dB)','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','bold');
legend({'N=64','N=128','N=256'},'FontName','Times New Roman','FontSize',12,'FontWeight','bold');
% filename = 'taylor.wav';
% audiowrite(filename,Reconstructed_array,Fs);
% filename = 'Dominque.wav';
% audiowrite(filename,Dom_recon_array,Fs);



