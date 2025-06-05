% Function to create Regression Matrix
function X = create_regression_matrix(y, u, nf, nb, nk)
    N = length(y);
    X = zeros(N, nf + nb); % Initialize Regression matrix

    for i = 1:N
        for j = 1:nf
            if (i - j) > 0
                X(i, j) = -y(i - j); 
            end
        end
        for j = 1:nb
            if (i - nk - j + 1) > 0
                X(i, nf + j) = u(i - nk - j + 1);
            end
        end
    end
end

% Function to estimate OE model parameters using Least Squares
function theta = estimate_oe(X, y)
    theta = (X' * X) \ (X' * y);
end

% Function to simulate OE model output


function y_oe = simulate_oe(u, theta, nf, nb, nk)
    disp(['Length of theta: ', num2str(length(theta))]);
    N = length(u);
    y_oe = zeros(N, 1);

    F_coeffs = [1; theta(1:nf)];
    B_coeffs = theta(nf+1:end);

    for i = 1:N
        for j = 1:min(nf, i-1)
        %for j = 1:nf
            if (i - j) > 0
                y_oe(i) = y_oe(i) - F_coeffs(j + 1) * y_oe(i - j);
            end
        end
        for j = 1:min(nb, i-nk)
        %for j = 1:nb
            if (i - nk - j + 1) > 0
                y_oe(i) = y_oe(i) + B_coeffs(j) * u(i - nk - j + 1);
            end
        end
    end
end

% Cost function for nonlinear optimization
function J = oe_cost_function(theta, y, u, nf, nb, nk)
    y_oe = simulate_oe(y, u, theta, nf, nb, nk);
    error = y - y_oe;
    J = sum(error.^2);
end

% Function to optimize OE model parameters using nonlinear optimization
function theta_opt = optimize_oe_parameters(y, u, nf, nb, nk, theta_init)
    options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter');
    theta_opt = fminunc(@(theta) oe_cost_function(theta, y, u, nf, nb, nk), theta_init, options);
end

% Define OE Model Orders
nf = 3;
nb = 3;
nk = 10;

% Load input-output signals
u = prbs_input;
y = y_prbs_G2;

% Least Squares Estimation
X_oe = create_regression_matrix(y, u, nf, nb, nk);
theta_oe = estimate_oe(X_oe, y);
y_oe = simulate_oe(u, theta_oe, nf, nb, nk);
RMSE_oe = sqrt(mean((y - y_oe).^2));
disp(['OE Model RMSE (Least Squares): ', num2str(RMSE_oe)]);

% Nonlinear Optimization
theta_init = rand(nf + nb, 1) * 0.001; % Small initial values
theta_opt = optimize_oe_parameters(y, u, nf, nb, nk, theta_init);
y_oe_opt = simulate_oe(y, u, theta_opt, nf, nb, nk);
RMSE_oe_opt = sqrt(mean((y - y_oe_opt).^2));
disp(['OE Model RMSE (Nonlinear Optimization): ', num2str(RMSE_oe_opt)]);

% Plot Comparison
time = 1:length(y);
figure;
plot(time, y, 'b', 'LineWidth', 1.5);
hold on;
plot(time, y_oe, 'r--', 'LineWidth', 1.5);
plot(time, y_oe_opt, 'g-.', 'LineWidth', 1.5);
hold off;
legend('True Output', 'OE Model (LS)', 'OE Model (Optimized)');
title('OE Model Comparison');
xlabel('Time Steps');
ylabel('Output');
grid on;