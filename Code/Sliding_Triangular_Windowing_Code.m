% Author: Chinmay Sahu
% Class: EE 501 - Digital Signal Processing
% Project: Signal Analysis/Synthesis using the DFT
% email: sahuc@clarkson.edu
% Clarkson Univeristy
% Nov 2017; Last revision: 19-Nov-2017
clc
clear
close all
filename = 'cleanspeech.wav';% Loading the sound file

Reconstructed_array=[]; % initializing reconstructed audio array for method-2
Dom_recon_array=[];  % initializing reconstructed audio array for method-1
first_n_picks=[]; %initializing array to keep first n picks for method-2
peak_keepers=[]; %initializing array to keep first n peaks for method-1


for fftPoints= 1 : 3
    
switch fftPoints
case 1
N=64;% Framerate
case 2
N=128;
case 3
N=256;
end


[audio_signal,Fs] = audioread('cleanspeech.wav');% Load the sound file


first_n_picks=zeros(1+(N/2),1);

for n=1:1:(1+(N/2)) % n varies from 1: 1+(N/2) example: for N=64, n=1:33
    
Reconstructed_array=[]; % clearing the reconstructed array for method-2 to start over again for each n 
Dom_recon_array=[]; % clearing the Dom_recon_array for method-2 to start over again for each n 

Snr_I_Num=0; % Snr_I_Num reset to zero to start over again for each n 
Snr_I_Den=0;  % Snr_I_Den reset to zero to start over again for each n 
Snr_II_Num=0;  % Snr_II_Num reset to zero to start over again for each n 
Snr_II_Den=0;  % Snr_II_Den reset to zero to start over again for each n 

Snr_Seg_I=0; % segmental SNR reset to zero to start over again for each n
Snr_Seg_II=0; % segmental SNR reset to zero to start over again for each n

 %
 %%
 % *Signal Synthesis:* For Sliding triangular windowing,audio signal
 % increments in (N/2)interval till {length of audio-(N/2)}
for li=1:N/2:(32768-(N/2)) % li: increments at N/2 intervals
    
Xt=audio_signal(li:li+(N-1)); % Signal that is selected from 1:64, 33:96, 65:127,
Y_1=Xt.*triang(N); % multiplying triangular window with each frame

N_frame_FFT=fft(Y_1); % N frame FFT of audio signal

% Method 1 dominant signal pick
[sortedValues,sortIndex] = sort(abs(N_frame_FFT(1:1+(N/2))),'descend');  %# Sort the values in descending order                                      
dominant_peaks = sortIndex(1:n); % n= number of dominant peaks
peak_keepers=zeros(1+(N/2),1); %initializing the peak_keepers to zeros 

for q=1:1:n
peak_keepers(dominant_peaks(q))=N_frame_FFT(dominant_peaks(q)); %
end

Dom_FlipSig=flipud(conj(peak_keepers));
Dom_Sym_Sig=[peak_keepers;Dom_FlipSig(2:N/2)];
Dom_N_frame_IFFT=ifft(Dom_Sym_Sig,N);


if li==1
Dom_recon_array=[Dom_recon_array;Dom_N_frame_IFFT];
elseif li ~=1
Pad_Dom_N_frame_IFFT=vertcat(zeros(li-1,1),Dom_N_frame_IFFT); %li=33,65,97
Dom_recon_array=vertcat(Dom_recon_array,zeros(N/2,1));
Dom_recon_array=Dom_recon_array+Pad_Dom_N_frame_IFFT;    
end


%Method-2 the first n components are retained 

first_n_picks=N_frame_FFT(1:n);
first_n_picks=[first_n_picks;zeros((1+(N/2)-n),1)];
FlipSig=flipud(conj(first_n_picks));
Sym_Sig=[first_n_picks;FlipSig(2:N/2)];
N_frame_IFFT=ifft(Sym_Sig,N);

if li==1
Reconstructed_array=[Reconstructed_array;N_frame_IFFT];
elseif li ~=1
Pad_N_frame_IFFT=vertcat(zeros(li-1,1),N_frame_IFFT); %li=33,65,97
Reconstructed_array=vertcat(Reconstructed_array,zeros(N/2,1));
Reconstructed_array=Reconstructed_array+Pad_N_frame_IFFT;    
end

end


Reconstructed_audio_array_I{n}=Dom_recon_array; % capturing audio signals for method-1
Reconstructed_audio_array_II{n}=Reconstructed_array;% capturing audio signals for method-2


%SNR_New
[SNR_I(n),SNR_SEGMENTAL_I(n),SNR_II(n),SNR_SEGMENTAL_II(n)] = CalculateSNR(audio_signal,Reconstructed_audio_array_I{n},Reconstructed_audio_array_II{n},N);

end % end of for loop where we move n using Method-1 or Method-2 to find SNR ,audios



switch fftPoints
case 1 % retrive audios and SNR for 64 frame
    
%Retrived Audio Signals
T_Compressed_Audio_Signals_64_1=Reconstructed_audio_array_I;
T_Compressed_Audio_Signals_64_2=Reconstructed_audio_array_II;   
% Overall SNR    
T_SNR_64_1=SNR_I;
T_SNR_64_2=SNR_II;
% Segmental
T_SNR_SEGMENTAL_64_1=SNR_SEGMENTAL_I;
T_SNR_SEGMENTAL_64_2=SNR_SEGMENTAL_II;

case 2
    
%Retrived Audio Signals
T_Compressed_Audio_Signals_128_1=Reconstructed_audio_array_I;
T_Compressed_Audio_Signals_128_2=Reconstructed_audio_array_II;  
% Overall SNR    
T_SNR_128_1=SNR_I;
T_SNR_128_2=SNR_II;
% Segmental
T_SNR_SEGMENTAL_128_1=SNR_SEGMENTAL_I;
T_SNR_SEGMENTAL_128_2=SNR_SEGMENTAL_II;

case 3
    
%Retrived Audio Signals
T_Compressed_Audio_Signals_256_1=Reconstructed_audio_array_I;
T_Compressed_Audio_Signals_256_2=Reconstructed_audio_array_II;      
% Overall SNR    
T_SNR_256_1=SNR_I;
T_SNR_256_2=SNR_II;
% Segmental
T_SNR_SEGMENTAL_256_1=SNR_SEGMENTAL_I;
T_SNR_SEGMENTAL_256_2=SNR_SEGMENTAL_II;
end

end


% Overall SNR method 1
figure
grid on
box on
grid minor
hold on
plot(100*((1:33)/64),T_SNR_64_1,':k','LineWidth',2);
hold on
plot(100*((1:65)/128),T_SNR_128_1,'--b','LineWidth',2);
hold on
plot(100*((1:129)/256),T_SNR_256_1,':r','LineWidth',2);
hold on
title('Sliding Triangular Method-1,Overall SNR vs % of components selected (n/N)','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
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
plot(100*((1:33)/64),T_SNR_SEGMENTAL_64_1,':k','LineWidth',2);
hold on
plot(100*((1:65)/128),T_SNR_SEGMENTAL_128_1,'--b','LineWidth',2);
hold on
plot(100*((1:129)/256),T_SNR_SEGMENTAL_256_1,':r','LineWidth',2);


title('Sliding Triangular Method-1,Segmental SNR vs % of components selected (n/N)','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
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
plot(100*((1:33)/64),T_SNR_64_2,':k','LineWidth',2);
hold on
plot(100*((1:65)/128),T_SNR_128_2,'--b','LineWidth',2);
hold on
plot(100*((1:129)/256),T_SNR_256_2,':r','LineWidth',2);


title('Sliding Triangular Method-2,Overall SNR vs % of components selected (n/N) ','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
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
plot(100*((1:33)/64),T_SNR_SEGMENTAL_64_2,':k','LineWidth',2);
hold on
plot(100*((1:65)/128),T_SNR_SEGMENTAL_128_2,'--b','LineWidth',2);
hold on
plot(100*((1:129)/256),T_SNR_SEGMENTAL_256_2,':r','LineWidth',2);

title('Sliding Triangular Method-2,Segmental SNR vs % of components selected (n/N)','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
xlabel('% of components selected (n/N)','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
ylabel('Segmental SNR(dB)','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','bold');
legend({'N=64','N=128','N=256'},'FontName','Times New Roman','FontSize',12,'FontWeight','bold');

% filename = 'taylor.wav';
% audiowrite(filename,Reconstructed_array,Fs);
% filename = 'Dominque.wav';
% audiowrite(filename,Dom_recon_array,Fs);



