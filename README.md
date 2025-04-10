Distance snr in awgn-
clc;clear all;close all;
pt = 1;                % Transmit power (W)
n = 2;                 % Path loss exponent
d = linspace(1, 10, 100);  % Distance from 1 to 10
% Calculate received signal power using path loss model
Pr = pt ./ (d.^n);     % Received signal power (deterministic)
% Plotting
plot(d, Pr, '-b', 'LineWidth', 2);
xlabel('Distance');
ylabel('Received Signal Power (W)');
title('Distance vs Received Signal Power');
grid on;

link capacity-
clc;clear;
snr_db = linspace(0, 20, 100);           % SNR in dB
snr = 10.^(snr_db / 10);                 % Convert to linear
capacity = log2(1 + snr);                % Shannon formula
plot(snr_db, capacity, 'r', 'LineWidth', 2);
xlabel('SNR (dB)');
ylabel('Channel Capacity (bps/Hz)');
title('AWGN Channel Capacity vs SNR');
grid on;

Okamura-
clc; clear;
% User inputs
fc = input('Enter frequency (MHz): ');
Gt = input('Enter transmitter antenna gain Gt (dB): ');
Gr = input('Enter receiver antenna gain Gr (dB): ');
kc = input('Enter correction factor kc (dB): ');
d_start = input('Enter starting distance (km): ');
d_end = input('Enter ending distance (km): ');
d = d_start:1:d_end; % distance vector
% Free space path loss
Lf = 32.45 + 20*log10(fc) + 20*log10(d);
% Median attenuation (simplified approximation)
amu = 10*log10(d); 
% Okumura Path Loss
for i = 1:length(d)
    l50freq(i) = Lf(i) + amu(i) - Gt - Gr - kc;
end
% Plotting
plot(d, l50freq, 'k-o');
xlabel('Distance (km)');
ylabel('Path Loss (dB)');
title('Okumura Model Path Loss with User Inputs');
grid on;

hata-
clc; clear;
% User inputs
fc = input('Enter frequency (MHz): '); % Frequency in MHz
hb = input('Enter base station height hb (meters): '); % Base station height in meters
hm = input('Enter mobile station height hm (meters): '); % Mobile station height in meters
d_start = input('Enter starting distance (km): '); % Starting distance
d_end = input('Enter ending distance (km): '); % Ending distance
d = d_start:1:d_end; % Distance vector
% Hata Model Path Loss for Urban Area
L_hata = 69.55 + 26.16*log10(fc) - 13.82*log10(hb) - (1.1*log10(fc) - 0.7)*hm + (44.9 - 6.55*log10(hb))*log10(d);
% Plotting
plot(d, L_hata, 'b-o');
xlabel('Distance (km)');
ylabel('Path Loss (dB)');
title('Hata Model Path Loss with User Inputs (Urban Area)');
grid on;

bpsk awgn-
clc; clear;
SNR = 0:10;
BER = zeros(size(SNR));
for i = 1:length(SNR)
    EbN0 = 10^(SNR(i)/10);
    noise_var = 0.5 / EbN0;
    bits = randi([0 1], 1, 1e6);
    s = 2*bits - 1;
    r = s + sqrt(noise_var) * randn(1, length(s));
    b_hat = (r < 0);
    BER(i) = sum(bits ~= b_hat) / length(bits);
end
semilogy(SNR, BER); grid on;
xlabel('E_b/N_0 (dB)'); ylabel('BER');
title('BPSK over AWGN');

qpsk awgn-
clc; clear;
SNR = 0:2:10;
BER = zeros(size(SNR));
for i = 1:length(SNR)
    EbN0 = 10^(SNR(i)/10);
    N = 1e6;
    sI = 2*randi([0 1],1,N)-1;
    sQ = 2*randi([0 1],1,N)-1;
    noise = (1/sqrt(2*EbN0))*(randn(1,N)+1j*randn(1,N));
    r = sI + 1j*sQ + noise;
    bI = real(r) < 0;
    bQ = imag(r) < 0;
    BER(i) = (sum(bI ~= (sI<0)) + sum(bQ ~= (sQ<0)))/(2*N);
end
semilogy(SNR, BER); grid on;
xlabel('E_b/N_0 (dB)'); ylabel('BER');
title('QPSK over AWGN');

bpsk Laplacian-
clc; clear;
SNR = 0:10;
BER = zeros(size(SNR));
for i = 1:length(SNR)
    EbN0 = 10^(SNR(i)/10);
    a = sqrt(0.25 / EbN0); 
    bits = randi([0 1], 1, 1e6);
    s = 2*bits - 1;
    u = rand(1, length(s)) - 0.5;
    n = a * sign(u) .* log(1 - 2*abs(u));
    r = s + n;
    b_hat = (r < 0);
    BER(i) = sum(bits ~= b_hat) / length(bits);
end
semilogy(SNR, BER); grid on;
xlabel('E_b/N_0 (dB)'); ylabel('BER');
title('BPSK over Laplacian');

qpsk Laplacian-
clc; clear;
SNR = 0:2:10;
BER = zeros(size(SNR));
for i = 1:length(SNR)
    EbN0 = 10^(SNR(i)/10);
    a = sqrt(0.25 / EbN0);
    N = 1e6;
    sI = 2*randi([0 1],1,N)-1;
    sQ = 2*randi([0 1],1,N)-1;
    uI = rand(1,N)-0.5;
    uQ = rand(1,N)-0.5;
    nI = a*sign(uI).*log(1-2*abs(uI));
    nQ = a*sign(uQ).*log(1-2*abs(uQ));
    r = (sI + nI) + 1j*(sQ + nQ);
    bI = real(r) < 0;
    bQ = imag(r) < 0;
    BER(i) = (sum(bI ~= (sI<0)) + sum(bQ ~= (sQ<0))) / (2*N);
end
semilogy(SNR, BER); grid on;
xlabel('E_b/N_0 (dB)'); ylabel('BER');
title('QPSK over Laplacian');
