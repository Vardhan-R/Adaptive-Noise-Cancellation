function metric = analyzeResults(x_mic, s_recovered, s_clean, fs, mode, tonalFreqs)
    N = length(x_mic);
    s_clean = s_clean(1:N);
    
    % SNR Evaluation (for full mode)
    if ~strcmp(mode, 'partial')
        noisy_SNR = 10*log10(sum(s_clean.^2) / sum((x_mic - s_clean).^2));
        enhanced_SNR = 10*log10(sum(s_clean.^2) / sum((s_recovered - s_clean).^2));
        metric = enhanced_SNR - noisy_SNR;  % SNR improvement
    else
        % Spectrum for partial mode (tone recovery)
        nfft = 8192;
        [Pxx_rec, f] = pwelch(s_recovered, hamming(1024), 512, nfft, fs);
        [Pxx_clean, ~] = pwelch(s_clean, hamming(1024), 512, nfft, fs);
        
        avg_tone_reduction = 0;
        for i = 1:length(tonalFreqs)
            center_freq = tonalFreqs(i);
            [~, idx] = min(abs(f - center_freq));
            tone_reduction_dB = 10*log10(Pxx_rec(idx)) - 10*log10(Pxx_clean(idx));
            avg_tone_reduction = avg_tone_reduction + tone_reduction_dB;
        end
        metric = avg_tone_reduction / length(tonalFreqs);  % Average tone recovery
    end
end
