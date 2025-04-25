% Load signals
% s_clean = load("clean_speech.txt");
% x_mic = load("noisy_speech.txt");
% w_ref = load("external_noise.txt");
[s_clean, fs] = audioread("audio_files/clnsp1.wav");
[x_mic, ~] = audioread("audio_files/noisy1_SNRdb_15.0_clnsp1.wav");
[w_ref, ~] = audioread("audio_files/noisy1_SNRdb_15.0.wav");

% Normalize signals
s_clean = s_clean / max(abs(s_clean));
x_mic = x_mic / max(abs(x_mic));
w_ref = w_ref / max(abs(w_ref));

% System parameters
% fs = 44100;
N = length(x_mic);
mode = 'full';
tonalFreqs = [1000, 2750, 5500, 8200];
r = 0.98;
filterOrder = 1;  % fixed

% === Golden Section Search Setup ===
a = 0.001;          % lower bound for mu
b = 1;              % upper bound for mu
tol = 1e-4;         % tolerance
gr = (sqrt(5) - 1) / 2;  % golden ratio coefficient

% Initial test points
c = b - gr * (b - a);
d = a + gr * (b - a);

% Evaluate function at initial points
snr_c = myLMS(fs, x_mic, w_ref, s_clean, tonalFreqs, c, filterOrder, mode);
snr_d = myLMS(fs, x_mic, w_ref, s_clean, tonalFreqs, d, filterOrder, mode);

% Iteratively narrow the interval
while (b - a) > tol
    if snr_c > snr_d
        b = d;
        d = c;
        snr_d = snr_c;
        c = b - gr * (b - a);
        snr_c = myLMS(fs, x_mic, w_ref, s_clean, tonalFreqs, c, filterOrder, mode);
    else
        a = c;
        c = d;
        snr_c = snr_d;
        d = a + gr * (b - a);
        snr_d = myLMS(fs, x_mic, w_ref, s_clean, tonalFreqs, d, filterOrder, mode);
    end
end

% Best mu is midpoint of final interval
best_mu = (a + b)/2;
best_snr = max(snr_c, snr_d);

% Output results
fprintf('Best μ (via GSS): %.6f\n', best_mu);
fprintf('Fixed Filter Order: %d\n', filterOrder);
fprintf('Maximum SNR: %.2f dB\n', best_snr);
