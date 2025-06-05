%% System Identification: OE Model (Simplified)
% Uses predefined systems (G1, G2), input signals, and parameters.
% ------------------------------------------------------------------------

%% Step 1: Estimate ARX Model for Initialization
% ------------------------------------------------------------------------
% Define model orders for System 2 (PT2-Tt)
nf = 2;     % Denominator order (F(q)) 
nb = 2;     % Numerator order (B(q)) 
nk = 10;    % Dead time in samples (1s / T0 = 10)

% Create valid ARX regression matrix
[X_arx, y_trimmed] = create_regression_matrix(y_prbs_G2, prbs_input, nf, nb, nk);

% Estimate ARX parameters
theta_arx = X_arx \ y_trimmed;

% Initialize OE parameters with ARX coefficients
theta_init_oe = [theta_arx(1:nf); theta_arx(nf+1:end)];

% Optimize OE model
theta_opt_oe = optimize_oe(y_prbs_G2, prbs_input, nf, nb, nk, theta_init_oe);

% Simulate OE output
y_oe = simulate_oe(prbs_input, theta_opt_oe, nf, nb, nk);
%% Step 2: Optimize OE Model with ARX Initialization
% ------------------------------------------------------------------------
% Extract ARX coefficients for OE initialization
theta_init_oe = [theta_arx(1:nf); theta_arx(nf+1:end)];

% Optimize OE parameters
theta_opt_oe = optimize_oe(y_prbs_G2, prbs_input, nf, nb, nk, theta_init_oe);

% Simulate OE model output
y_oe = simulate_oe(prbs_input, theta_opt_oe, nf, nb, nk);

%% Step 3: Validation and Plots
% ------------------------------------------------------------------------
% Calculate RMSE (compare with noiseless output)
y_true = lsim(G2, prbs_input, t); % True (noiseless) output
RMSE_oe = sqrt(mean((y_true - y_oe).^2));

% Plot results
figure;
plot(t, y_prbs_G2, 'b', 'LineWidth', 1.5); hold on;
plot(t, y_oe, 'r--', 'LineWidth', 1.5); hold off;
xlim([0, 20]); % Zoom into first 20 seconds
legend('Measured Output (Noisy)', 'OE Model');
title('OE Model Validation (System 2)');
xlabel('Time (s)'); ylabel('Output'); grid on;

% Display results
disp('=======================');
disp('OE Model Results');
disp('=======================');
disp(['RMSE: ', num2str(RMSE_oe)]);
disp('Estimated F(q) = 1 + '); disp(theta_opt_oe(1:nf)');
disp('Estimated B(q) = '); disp(theta_opt_oe(nf+1:end)');

%% Supporting Functions
% ------------------------------------------------------------------------
function [X, y_trimmed] = create_regression_matrix(y, u, na, nb, nk)
    % Ensure column vectors and valid indices
    y = y(:);
    u = u(:);
    N = length(y);
    
    % Determine valid start and end indices
    start_idx = max(na, nb + nk) + 1;
    end_idx = N;
    num_rows = end_idx - start_idx + 1;
    
    % Initialize regression matrix
    X = zeros(num_rows, na + nb);
    
    % Fill matrix with valid indices only
    for i = 1:num_rows
        current_idx = start_idx + i - 1;
        
        % Past outputs (AR part)
        X(i, 1:na) = -y(current_idx - 1 : -1 : current_idx - na);
        
        % Past inputs (X part)
        X(i, na+1:end) = u(current_idx - nk : -1 : current_idx - nk - nb + 1);
    end
    
    % Trim output to match regression matrix
    y_trimmed = y(start_idx:end);
end

function y_oe = simulate_oe(u, theta, nf, nb, nk)
    N = length(u);
    y_oe = zeros(N, 1);
    F = [1; theta(1:nf)];
    B = theta(nf+1:end);
    
    for i = 1:N
        % Feedback terms (F(q))
        for j = 1:min(nf, i-1)
            y_oe(i) = y_oe(i) - F(j+1) * y_oe(i-j);
        end
        
        % Feedforward terms (B(q))
        for j = 1:min(nb, i-nk)
            if (i - nk - j + 1) > 0
                y_oe(i) = y_oe(i) + B(j) * u(i - nk - j + 1);
            end
        end
    end
end

function theta_opt = optimize_oe(y, u, nf, nb, nk, theta_init)
    function cost = oe_cost(theta)
        y_oe = simulate_oe(u, theta, nf, nb, nk);
        valid_indices = max(nf, nb + nk) + 1 : length(y);
        cost = sum((y(valid_indices) - y_oe(valid_indices)).^2);
    end
    options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter');
    theta_opt = fminunc(@oe_cost, theta_init, options);
end