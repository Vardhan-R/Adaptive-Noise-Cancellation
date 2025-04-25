% Load signals
s_clean = load("clean_speech.txt");
x_mic = load("noisy_speech.txt");
w_ref = load("external_noise.txt");

% % Normalize signals
% x_mic = x_mic / max(abs(x_mic));
% w_ref = w_ref / max(abs(w_ref));

% System parameters
fs = 44100;
N = length(x_mic);
mode = 'partial';
tonalFreqs = [1000, 2750, 5500, 8200];
r = 0.99;
filterOrder = 5;  % Fix filter order for now

% === Golden Section Search for lambda ===
a = 0.1;
b = 0.9999;
tol = 1e-3;
gr = (sqrt(5) - 1) / 2;

% Initial points
c = b - gr * (b - a);
d = a + gr * (b - a);

% Evaluate at initial points
% snr_c = myRLS(x_mic, w_ref, s_clean, tonalFreqs, c, filterOrder, mode);
% snr_d = myRLS(x_mic, w_ref, s_clean, tonalFreqs, d, filterOrder, mode);
[s_rec_c, ~] = adaptiveNoiseFilter(x_mic, w_ref, filterOrder, fs, mode, tonalFreqs, r, c);
[s_rec_d, ~] = adaptiveNoiseFilter(x_mic, w_ref, filterOrder, fs, mode, tonalFreqs, r, d);
snr_c = analyzeResults(x_mic, s_rec_c, s_clean, fs, mode, tonalFreqs);
snr_d = analyzeResults(x_mic, s_rec_d, s_clean, fs, mode, tonalFreqs);

% Iterative narrowing
while (b - a) > tol
    if snr_c > snr_d
        b = d;
        d = c;
        snr_d = snr_c;
        c = b - gr * (b - a);
        % snr_c = myRLS(x_mic, w_ref, s_clean, tonalFreqs, c, filterOrder, mode);
        [s_rec_c, ~] = adaptiveNoiseFilter(x_mic, w_ref, filterOrder, fs, mode, tonalFreqs, r, c);
        snr_c = analyzeResults(x_mic, s_rec_c, s_clean, fs, mode, tonalFreqs);
    else
        a = c;
        c = d;
        snr_c = snr_d;
        d = a + gr * (b - a);
        % snr_d = myRLS(x_mic, w_ref, s_clean, tonalFreqs, d, filterOrder, mode);
        [s_rec_d, ~] = adaptiveNoiseFilter(x_mic, w_ref, filterOrder, fs, mode, tonalFreqs, r, d);
        snr_d = analyzeResults(x_mic, s_rec_d, s_clean, fs, mode, tonalFreqs);
    end
    
end

% Best lambda
best_lambda = (a + b) / 2;
best_snr = max(snr_c, snr_d);

% Display
fprintf('Mode: %s\n', mode);
fprintf('Fixed r: %d\n', r);
fprintf('Fixed Filter Order: %d\n', filterOrder);
fprintf('Best λ (via GSS): %.6f\n', best_lambda);
if strcmp(mode, 'partial')
    fprintf('Average tone recovery: %.2f dB\n', best_snr);
else
    fprintf('Maximum SNR: %.2f dB\n', best_snr);
end
