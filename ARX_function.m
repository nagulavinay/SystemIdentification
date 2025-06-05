% Function to create Regression Matrix
function X = create_regression_matrix(y, u, na, nb, nk)
        % y: Output data
        % u: Input data
        % na: Order of A(q)
        % nb: Order of B(q)
        % nk: Dead time (delay)
 N = length(y); % No. of Samples
 X = zeros(N, na + nb); % Initialize Regression matrix

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
end

% Function to compute model output
function y_model = simulate_arx(y, u, theta, na, nb, nk)
    N = length(y);
    y_model = zeros(N, 1);
    %y_model = zeros(size(y));
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
rng(123);
na = 3;
nb = 1;
nk = 0;
u = step_input
y = y_step_G1
X = create_regression_matrix(y, u, na, nb, nk);
theta = estimate_arx (X, y);
y_arx = simulate_arx(y, u, theta, na, nb, nk);
disp('Before rec filtering');
disp(y_arx);
figure;
plot(t, y, 'b', t, y_arx, 'r--');
legend('True Output', 'ARX Model Output');
title('ARX Model vs. True Output');
xlabel('Time(s)');
ylabel('Output');
grid on;

% Compute RMSE for Model validation
RMSE = sqrt(mean((y - y_arx).^2));
disp(['Step Input (G1) - ARX Model RMSE: ', num2str(RMSE)]);
figure;
plot(1:length(y), y, 'b', 1:length(y_arx), y_arx, 'r--');
legend('True Output', 'ARX Model Output');
title('ARX Model vs. True Output - Step Input (G1)');
xlabel('Time steps');
ylabel('Output');
grid on;

% Instrumental Variable Method
function theta_iv = estimate_arx_iv(y, u, z, na, nb, nk)
    X = create_regression_matrix(y, u, na, nb, nk);
    Z = create_regression_matrix(y, z, na, nb, nk);
    theta_iv = (Z' * X) \ (Z' * y);
end
z = filter([0.5 -0.5], [1 -0.8], u);
z = [zeros(nk, 1); z(1:end-nk)];
theta_iv = estimate_arx_iv(y, u, z, na, nb, nk);
y_arx_iv = simulate_arx(y, u, theta_iv, na, nb, nk);
RMSE_iv = sqrt(mean((y - y_arx_iv).^2));
disp(['Instrumental Variables ARX Model RMSE: ', num2str(RMSE_iv)]);
figure;
plot(1:length(y), y, 'b', 1:length(y_arx_iv), y_arx_iv, 'r--');
legend('True Output', 'ARX Model Output with IV');
title('ARX Model with Instrumental Variables');
xlabel('Time Steps');
ylabel('Output');
grid on;

% Estimation using Recursive filtering
X = create_regression_matrix(y, u, na, nb, nk);
theta = estimate_arx(X, y);

% Step 2: Extract denominator coefficients (F polynomial)
F = [1; theta(1:na).']; % Ensure column vector
%F = F / F(1); % Normalizing first coeffecient
if any(abs(roots(F)) >= 1)
    error('Unstable filter detected! Recursive filtering may be inaccurate.');
end
% Step 3: Filter the input and output using 1/F(q)
y_filt = filter(1, F, y);  % Apply 1/F filtering to output
u_filt = filter(1, F, u);  % Apply 1/F filtering to input

% Step 4: Re-estimate ARX model with filtered data
X_filt = create_regression_matrix(y_filt, u_filt, na, nb, nk);
theta_filt = estimate_arx(X_filt, y_filt);

% Step 5: Simulate the filtered ARX model output
y_arx_filt = simulate_arx(y_filt, u_filt, theta_filt, na, nb, nk);
disp('After Recursive filtering');
disp(y_arx_filt)



% Plot comparison
N = min(length(y), length(y_arx_filt)); % Find shortest length
time = 1:N; % Create time vector based on the common length
%time_filtered = time(2:end); % Exclude the first value for filtering
%N = min(length(y_filt), length(y_arx_filt)); % Shortest length
time_filtered = 1:N;
y_filtered = y(1:N);
y_arx_filt_filtered = y_arx_filt(1:N);

% Compute RMSE for Recursive Filtering ARX Model
RMSE_filt = sqrt(mean((y_filt(1:N) - y_arx_filt(1:N)).^2));
disp(['Recursive Filtering ARX Model RMSE: ', num2str(RMSE_filt)]);
% Plot comparison (excluding first value)
figure;
plot(time_filtered, y_filtered, 'b', 'LineWidth', 1.5); % True Output
hold on;
plot(time_filtered, y_arx_filt_filtered, 'r--', 'LineWidth', 1.5); % Filtered ARX Output
hold off;
legend('True Output', 'ARX Model Output with Filtering');
title('ARX Model with Recursive Filtering');
xlabel('Time Steps');
ylabel('Output');
grid on;

