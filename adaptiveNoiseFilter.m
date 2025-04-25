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
