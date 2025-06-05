function X = create_fir_regression_matrix(u, nb_fir, nk)
        % Function to construct the FIR regression matrix
        % u: Input signal
        % nb: Order of FIR model
        % nk: Dead time (delay)
    N = length(u);
    X = zeros(N, nb_fir);
    % Construct FIR regression matrix using past input values
    for i = 1:N
        for j = 1:nb_fir
            if (i - nk - j + 1) > 0
                X(i, j) = u(i - nk - j + 1);
            end
        end
    end
end

% Function to estimate FIR coeffecients using LS
function theta_fir = estimate_fir(X, y)

    theta_fir = (X' * X) \ (X' * y);
end

% Function to estimate FIR model's predicted output
function y_fir = simulate_fir(u, theta_fir, nb, nk)
        % u: Input signal
        % theta_fir: FIR coefficients
        % nb: Order of FIR model
        % nk: Dead time
    N = length(u);
    y_fir = zeros(N, 1);
    for i = 1:N
        for j = 1:nb
            if (i - nk - j + 1) > 0
                y_fir(i) = y_fir(i) + theta_fir(j) * u(i - nk - j +1);
            end
        end
    end
end

u = step_input;
y = y_step_G1;
nb_fir = 30;
nk = 0;
X_fir = create_fir_regression_matrix(u, nb_fir, nk);
theta_fir = estimate_fir(X_fir, y);
y_fir = simulate_fir(u, theta_fir, nb_fir, nk);
N_min3 = min(length(y), length(y_fir));
y_fil = y(1:N_min3);
y_fir = y_fir(1:N_min3);
RMSE_fir = sqrt(mean((y_fil(:) - y_fir(:)).^2));
disp(['Step Input (G1) - FIR Model RMSE: ', num2str(RMSE_fir)]);

figure;
%h{1} = plot(1:length(y), y, 'b', 'LineWidth', 1.5);
plot(1:N_min3, y_fil, 'b', 'LineWidth', 1.5);
hold on;
%h{2} = plot(1:length(y_fir), y_fir, 'r--', 'LineWidth', 1.5);
plot(1:N_min3, y_fir, 'r--', 'LineWidth', 1.5);

hold off;
legend('True Output', 'FIR Model Output', 'Location', 'northeast');
title('FIR estimation for step signal');
xlabel('K samples');
ylabel('Output');
grid on;