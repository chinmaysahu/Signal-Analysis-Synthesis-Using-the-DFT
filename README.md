# Signal-Analysis-Synthesis-Using-the-DFT

This project presents a discrete Fourier transform (DFT) based
audio compression technique.This involves a three-step process.
First, the audio signal is transformed in terms of frequency
components by using DFT. Next, the audio signal
is compressed by selectively picking peaks in DFT magnitude
spectrum while removing redundant elements. Peak
picking is performed either by selecting n dominant components(
Method 1) or by choosing first n components (Method
2). Finally, the resultant frequency domain results are transformed
back to time domain using inverse DFT. Sliding
rectangular window and sliding triangular window are multiplied
to the audio signal for different values of N-point FFT
to evaluate the quality of reconstructed audio. Using reconstructed
audio signal, overall and segmental signal to noise
ratio(SNR) is estimated for different values of N(64,128,256).
The results are tabulated, analyzed based on different values
of peak picking(n).
