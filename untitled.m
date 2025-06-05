% Function to create Regression Matrix
function X = create_regression_matrix(y, u, na, nb, nk)
        % y: Output data
        % u: Input data
        % na: Order of A(q)
        % nb: Order of B(q)
        % nk: Dead time (delay)
 N = length(y); % No. of Samples
 X = zeros(N - na, na + nb); % Initialize Regression matrix

% Construct ARX regression matrix
    for i = 1:N
        % fill in past output values (AR part)
        for j = 1:na
            if (i - j) > 0
                X(i, j) = -y(i - j);
            end
        end
        % fill in past input values (X part)
        for j = 1:nb
            if(i -nk -j +1) > 0
                X(i, na + j) = u(i -nk - j + 1);
            end
        end
    end
end
% Function to estimate parameters using LS
function theta = estimate_arx(X,y)
         % X: Regression matrix
         % y: Output data
    theta = (X' * X) \ (X' * y);
    %disp(theta);
end

% Function to compute model output
function y_model = simulate_arx(y, u, theta, na, nb, nk)
    N = length(y);
    y_model = zeros(size(y));
    % compute simulated output recursively
    for i = 1:N
        for j = 1:na
            if (i - j) >0
                y_model(i) = y_model(i) - theta(j) * y_model(i - j);
            end
        end
        for j = 1:nb
            if (i - nk - j + 1) > 0
                y_model(i) = y_model(i) + theta(na + j) * u(i - nk - j +1);
            end
        end
    end
end

% Define ARX model orders
%rng(123);
na = 2;
nb = 1;
nk = 10;
u = prbs_input;
y = y_prbs_G2;
X = create_regression_matrix(y, u, na, nb, nk);
theta = estimate_arx (X, y);
y_arx = simulate_arx(y, u, theta, na, nb, nk);
%figure;
%plot(t, y, 'b', t, y_arx, 'r--');
%legend('True Output', 'ARX Model Output');
%title('ARX Model vs. True Output');
%xlabel('Time(s)');
%ylabel('Output');
%grid on;

% Compute RMSE for Model validation
N_min1 = min(length(y), length(y_arx));
y_f = y(1:N_min1);
y_arx = y_arx(1:N_min1);
RMSE = sqrt(mean((y_f(:) - y_arx(:)).^2));
disp(['Step Input (G1) - ARX Model RMSE: ', num2str(RMSE)]);

figure;
%h{1} = plot(1:length(y), y, 'b')
plot(1:N_min1, y_f, 'b', 'LineWidth', 1.5);
hold on;
plot(1:N_min1, y_arx, 'r--', 'LineWidth', 1.5);
%h{2} = plot(1:length(y_arx), y_arx, 'r--');
hold off;
%legend([h{1}(1);h{2}(1)], 'True Output', 'ARX Model Output');
title('ARX estimation (System 2)');
xlabel('K samples');
ylabel('Output');
legend('Actual output', 'Estimated ARX');
%grid on;



% Step 1: Extract F(q) = [1 + a1*q^-1 + ... + ana*q^-na]
F = [theta(1:na).'];
if any(abs(roots(F)) >= 1)
    error('Unstable recursive filter! Cannot apply 1/F(q).');
end
% Step 2: Apply 1/F(q) to y and u using filter
y_filt = filter(1, F, y);
u_filt = filter(1, F, u);

% Step 3: Create regression matrix with filtered signals
X_filt = create_regression_matrix(y_filt, u_filt, na, nb, nk);

% Step 4: Estimate new parameters from filtered signals
theta_filt = estimate_arx(X_filt, y_filt);

% Step 5: Simulate output using filtered model
y_arx_filt = simulate_arx(y_filt, u_filt, theta_filt, na, nb, nk);

% Make lengths consistent (truncate if needed)
N_min = min(length(y_filt), length(y_arx_filt));
y_filt_trimmed = y_filt(1:N_min);
y_arx_filt_trimmed = y_arx_filt(1:N_min);

% Step 6: Compute RMSE
RMSE_filt = sqrt(mean((y_filt_trimmed(:) - y_arx_filt_trimmed(:)).^2));
%fprintf('RMSE: %.4f\n', RMSE_iv);
disp(['Recursive Filtering ARX Model RMSE: ', num2str(RMSE_filt)]);

% Step 7: Plot
figure;
plot(1:N_min, y_filt_trimmed, 'b', 'LineWidth', 1.5);
hold on;
plot(1:N_min, y_arx_filt_trimmed, 'r--', 'LineWidth', 1.5);
hold off;
legend('True Output', 'Filtered ARX Output');
title('ARX Model with Recursive Filtering for step signal');
xlabel('K samples');
ylabel('Output');
grid off;






% Function to estimate using IV method
function Z = construct_instrumental_variables(y_arx, u, na, nb, nk)
    N = length(y_arx);
    Z = zeros(N, na + nb);
    for i = 1:N
        for j = 1:nb
            if (i - nk - j + 1) > 0
                Z(i, j) = u(i -nk - j + 1);
            end
        end
        for j = 1:na
            if(i - j) > 0
                Z(i, nb + j) = -y_arx(i - j);
            end
        end
    end
end

function theta_iv = estimate_arx_iv(y, u, na, nb, nk)
    X = create_regression_matrix(y,u, na, nb, nk);
    theta_arx = estimate_arx(X, y);
    num_iter = 2;
    for iter = 1:num_iter
        y_arx = simulate_arx(y, u, theta_arx, na, nb, nk);
        Z = construct_instrumental_variables(y_arx, u , na, nb, nk);
        theta_arx = (Z' * X) \ (Z' * y);
    end
    theta_iv = theta_arx;
end
theta_iv = estimate_arx_iv(y, u, na, nb, nk);
y_arx_iv = simulate_arx(y, u, theta_iv, na, nb, nk);
N_min2 = min(length(y), length(y_arx_iv));
y_ff = y(1:N_min);
y_arx_iv_f = y_arx_iv(1:N_min);
RMSE_iv = sqrt(mean((y_ff(:) - y_arx_iv_f(:)).^2));
%fprintf('RMSE: %.4f\n', RMSE_iv);
disp(['Instrumental Variables ARX Model RMSE: ', num2str(RMSE_iv)]);
figure;
%h{1} = plot(1:length(y), y, 'b');
plot(1:N_min, y_ff, 'b', 'LineWidth', 1.5);
hold on;
plot(1:N_min, y_arx_iv_f, 'r--', 'LineWidth', 1.5);
%h{2} = plot(1:length(y_arx), y_arx, 'r--');
%h{3} = plot(1:length(y_arx_iv), y_arx_iv, 'g-.');
hold off;
legend('True Output', 'IV Optimized ARX');
title('IV Optimized ARX Model for step signal');
xlabel('K samples');
ylabel('Output');
grid off;