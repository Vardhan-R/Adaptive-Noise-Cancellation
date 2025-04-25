% Main script for noise cancellation system
clc; clear; close all;

% Load signals
s_clean = load("clean_speech.txt");
x_mic = load("noisy_speech.txt");
w_ref = load("external_noise.txt");

% Parameters
fs = 44100;
filterOrder = 5;
mode = 'partial'; % 'full' or 'partial'
tonalFreqs = [1480, 2725, 5458]; % Hz
r = 0.999;
lambda = 0.99999;

% Perform adaptive filtering
[s_recovered, w_leak_hat] = adaptiveNoiseFilter(x_mic, w_ref, filterOrder, fs, mode, tonalFreqs, r, lambda);

% Save result
audiowrite('recovered_speech_rls.wav', s_recovered / max(abs(s_recovered)) * 0.9, fs);

% Analyze and visualize results
analyzeResults(x_mic, w_ref, s_recovered, w_leak_hat, s_clean, fs, mode, tonalFreqs);
% processAllTestCases();

%% Adaptive Noise Filtering Function
function [s_recovered, w_leak_hat] = adaptiveNoiseFilter(x_mic, w_ref, filterOrder, fs, mode, tonalFreqs, r, lambda)
    N = length(x_mic);
    
    % Initialization
    w_leak_hat = zeros(N, 1);
    s_recovered = zeros(N, 1);
    w_rls = zeros(filterOrder, 1);
    P = 1000 * eye(filterOrder);
    w_ref_buf = zeros(filterOrder, 1);
    
    % Notch filter setup (for partial mode)
    notch = {};
    if strcmp(mode, 'partial')
        for k = 1:length(tonalFreqs)
            omega0 = 2 * pi * tonalFreqs(k) / fs;
            b = [1, -2*cos(omega0), 1];
            a = [1, -2*r*cos(omega0), r^2];
            notch{k} = struct('b', b, 'a', a, 'z', zeros(2, 1));
        end
    end
    
    % Adaptive filtering
    for n = 1:N
        w_ref_n = w_ref(n);
        
        % Apply notch filters in partial mode
        if strcmp(mode, 'partial')
            for k = 1:length(notch)
                b = notch{k}.b; a = notch{k}.a; z = notch{k}.z;
                y_n = b(1)*w_ref_n + z(1);
                z(1) = b(2)*w_ref_n - a(2)*y_n + z(2);
                z(2) = b(3)*w_ref_n - a(3)*y_n;
                w_ref_n = y_n;
                notch{k}.z = z;
            end
        end
        
        % Update buffer
        w_ref_buf = [w_ref_n; w_ref_buf(1:end-1)];
        
        % RLS algorithm
        if n >= filterOrder
            x_vec = w_ref_buf;
            g = (P * x_vec) / (lambda + x_vec' * P * x_vec);
            w_leak_hat(n) = w_rls' * x_vec;
            e = x_mic(n) - w_leak_hat(n); % s + v - vhat = speaker output
            w_rls = w_rls + g * e;
            P = (P - g * x_vec' * P) / lambda;
            s_recovered(n) = e;
        else
            w_leak_hat(n) = 0;
            s_recovered(n) = x_mic(n);
        end
    end
end

%% Analysis and Visualization Function
function analyzeResults(x_mic, w_ref, s_recovered, w_leak_hat, s_clean, fs, mode, tonalFreqs)
    N = length(x_mic);
    s_clean = s_clean(1:N);
    
    % Print mode
    fprintf('Mode: %s\n', upper(mode));
    
    % SNR Evaluation (only for full mode)
    if ~strcmp(mode, 'partial')
        noisy_SNR = 10*log10(sum(s_clean.^2) / sum((x_mic - s_clean).^2));
        enhanced_SNR = 10*log10(sum(s_clean.^2) / sum((s_recovered - s_clean).^2));
        
        % Print SNR metrics
        fprintf('SNR Before Processing: %.2f dB\n', noisy_SNR);
        fprintf('SNR After Processing: %.2f dB\n', enhanced_SNR);
        fprintf('SNR Improvement: %.2f dB\n', enhanced_SNR - noisy_SNR);
    end
    
    % Spectrum analysis
    nfft = 8192;
    [Pxx_mic, f] = pwelch(x_mic, hamming(1024), 512, nfft, fs);
    [Pxx_ref, ~] = pwelch(w_ref, hamming(1024), 512, nfft, fs);
    [Pxx_rec, ~] = pwelch(s_recovered, hamming(1024), 512, nfft, fs);
    [Pxx_noise, ~] = pwelch(w_leak_hat, hamming(1024), 512, nfft, fs);
    [Pxx_clean, ~] = pwelch(s_clean, hamming(1024), 512, nfft, fs);
    
    % Plot spectrum
    figure;
    plot(f, 10*log10(Pxx_mic), 'k', 'LineWidth', 1.5); hold on;
    plot(f, 10*log10(Pxx_ref), 'm', 'LineWidth', 1.2);
    plot(f, 10*log10(Pxx_rec), 'b', 'LineWidth', 1.5);
    plot(f, 10*log10(Pxx_noise), 'r', 'LineWidth', 1.2);
    plot(f, 10*log10(Pxx_clean), 'g', 'LineWidth', 1.5);
    xlabel('Frequency (Hz)');
    ylabel('Power (dB)');
    title(['Power Spectrum - ' upper(mode) ' Mode']);
    legend({'Mic (Noisy)', 'Reference Noise', 'Recovered Speech', 'Estimated Noise', 'Clean Speech'});
    grid on;
    
    % Calculate tone reduction for partial mode
    if strcmp(mode, 'partial')
        fprintf('\nTone Reduction Analysis:\n');
        fprintf('----------------------\n');
        
        % Individual tone reduction
        avg_tone_reduction = 0;
        for i = 1:length(tonalFreqs)
            center_freq = tonalFreqs(i);
            [~, idx] = min(abs(f - center_freq));
            tone_reduction_dB = 10*log10(Pxx_rec(idx)) - 10*log10(Pxx_clean(idx));
            fprintf('Tone reduction at %d Hz: %.2f dB\n', center_freq, tone_reduction_dB);
            avg_tone_reduction = avg_tone_reduction + tone_reduction_dB;
        end
        
        % Average tone reduction
        avg_tone_reduction = avg_tone_reduction / length(tonalFreqs);
        fprintf('Average tone recovery across all frequencies: %.2f dB\n', avg_tone_reduction);
        
        % Add vertical lines at tonal frequencies on the spectrum plot
        figure(gcf); % Get current figure (the spectrum plot)
        for i = 1:length(tonalFreqs)
            xline(tonalFreqs(i), '--', ['f = ' num2str(tonalFreqs(i)) ' Hz'], 'LineWidth', 1.5);
        end
    end
end

function processAllTestCases()
    % Get list of test case folders
    baseFolder = 'test_cases';
    testFolders = dir(fullfile(baseFolder, 'case_*'));
    
    % Create results folder
    resultsFolder = fullfile(baseFolder, 'results');
    if ~exist(resultsFolder, 'dir')
        mkdir(resultsFolder);
    end
    
    % Process each test case
    fprintf('Processing test cases...\n');
    
    % Parameters for all test cases
    fs = 44100;
    filterOrder = 5;
    modes = {'full', 'partial'};
    tonalFreqs = [1000, 2725]; % Hz
    r = 0.999;
    lambda = 0.99999;
    
    % Create summary file
    summaryFile = fopen(fullfile(resultsFolder, 'summary_results.txt'), 'w');
    fprintf(summaryFile, 'Test Case Summary Results\n');
    fprintf(summaryFile, '======================\n\n');
    
    % Process each test case
    for i = 1:length(testFolders)
        testFolder = fullfile(baseFolder, testFolders(i).name);
        fprintf('Processing %s...\n', testFolders(i).name);
        
        % Load test data
        s_clean = load(fullfile(testFolder, 'clean_speech.txt'));
        x_mic = load(fullfile(testFolder, 'noisy_speech.txt'));
        w_ref = load(fullfile(testFolder, 'external_noise.txt'));
        
        % Load filter info
        filterInfo = fileread(fullfile(testFolder, 'filter_info.txt'));
        
        % Create result folder for this test case
        caseResultFolder = fullfile(resultsFolder, testFolders(i).name);
        if ~exist(caseResultFolder, 'dir')
            mkdir(caseResultFolder);
        end
        
        % Write original filter info to results
        fid = fopen(fullfile(caseResultFolder, 'filter_info.txt'), 'w');
        fprintf(fid, '%s\n', filterInfo);
        fclose(fid);
        
        % Process both modes
        for j = 1:length(modes)
            mode = modes{j};
            
            % Create mode-specific folder
            modeFolder = fullfile(caseResultFolder, mode);
            if ~exist(modeFolder, 'dir')
                mkdir(modeFolder);
            end
            
            % Process with adaptive filter
            [s_recovered, w_leak_hat] = adaptiveNoiseFilter(x_mic, w_ref, filterOrder, fs, mode, tonalFreqs, r, lambda);
            
            % Save audio result
            audiowrite(fullfile(modeFolder, 'recovered_speech.wav'), s_recovered / max(abs(s_recovered)) * 0.9, fs);
            
            % Save numeric result
            writematrix(s_recovered, fullfile(modeFolder, 'recovered_speech.txt'));
            
            % Redirect output to capture analysis results
            diary(fullfile(modeFolder, 'analysis_results.txt'));
            diary on;
            
            % Create figure for this case/mode
            figure;
            
            % Run analysis
            analyzeResults(x_mic, w_ref, s_recovered, w_leak_hat, s_clean, fs, mode, tonalFreqs);
            
            % Turn off output capture
            diary off;
            
            % Save figure
            saveas(gcf, fullfile(modeFolder, 'spectrum_analysis.png'));
            saveas(gcf, fullfile(modeFolder, 'spectrum_analysis.fig'));
            close(gcf);
            
            % Copy analysis results to summary
            analysisResults = fileread(fullfile(modeFolder, 'analysis_results.txt'));
            fprintf(summaryFile, '--- %s: %s MODE ---\n', testFolders(i).name, upper(mode));
            fprintf(summaryFile, '%s\n\n', analysisResults);
        end
    end
    
    fclose(summaryFile);
    fprintf('All test cases processed. Results saved in %s\n', resultsFolder);
end
