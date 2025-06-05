% Function to simulate OE output
function y_oe = simulate_oe(u, theta, nf, nb, nk)
    N = length(u);
    y_oe = zeros(N, 1);
    F = [1; theta(1:nf)];   % F(q) coefficients
    B = theta(nf+1:end);     % B(q) coefficients
    
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

% Function to optimize OE parameters
function theta_opt = optimize_oe(y, u, nf, nb, nk, theta_init)
    % Nested cost function
    function cost = oe_cost(theta)
        y_oe = simulate_oe(u, theta, nf, nb, nk);
        cost = sum((y - y_oe).^2);  % Sum of squared errors
    end
    
    % Optimize
    options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter');
    theta_opt = fminunc(@oe_cost, theta_init, options);
end

% Example usage
nf = 3;     % Denominator order (matches System 2: s^2 + 1.5s + 0.5)
nb = 1;     % Numerator order (matches System 2: 0.25)
nk = 0;    % Dead time in samples (e.g., T0 = 0.1s, dead time = 1s → nk=10)

% Load/generate input-output data
u = step_input(:);  % PRBS input (column vector)
y = y_step_G1(:);   % True output (column vector)

% Ensure u and y have the same length
min_length = min(length(u), length(y));
u = u(1:min_length);
y = y(1:min_length);

% Initial guess for parameters
theta_init = rand(nf + nb, 1); 

% Optimize
theta_opt = optimize_oe(y, u, nf, nb, nk, theta_init);

% Simulate and plot
y_oe = simulate_oe(u, theta_opt, nf, nb, nk);
RMSE = sqrt(mean((y - y_oe).^2));
disp(['RMSE: ', num2str(RMSE)]);

figure;
plot(y, 'b', 'LineWidth', 1.5); hold on;
plot(y_oe, 'r--', 'LineWidth', 1.5); hold off;
legend('True Output', 'OE Model');
title('OE Model optimization for step response');
xlabel('K samples'); ylabel('Output'); grid off;