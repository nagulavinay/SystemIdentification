X = create_regression_matrix(y, u, na, nb, nk);
theta_arx = estimate_arx(X, y);
B = [1; theta_arx(1:na)];
F = theta_arx(na+1:end);
y_filt = filter(1, F, y); 
u_filt = filter(1, F, u);
X_filt = create_regression_matrix(y_filt, u_filt, na, nb, nk);
theta_filt = estimate_arx(X_filt, y_filt);
F_oe = [1; theta_filt(1:na)];
B_oe = theta_filt(na+1:end);
y_oe = simulate_oe(u, [F_oe(2:end); B_oe], na, nb, nk);
function y_oe = simulate_oe(u, theta, nf, nb, nk)
    N = length(u);
    y_oe = zeros(N,1);
    B = [1; theta(1:nf)];
    F = theta(nf+1:end);
    
    for k = 1:N
        for i = 1:nf
            if k - i > 0
                y_oe(k) = y_oe(k) - F(i+1) * y_oe(k - i);
            end
        end
        for i = 1:nb
            if k - nk - i + 1 > 0
                y_oe(k) = y_oe(k) + B(i) * u(k - nk - i + 1);
            end
        end
    end
end
na = 2;
nb = 1; nf = 2;
RMSE_oe = sqrt(mean((y - y_oe).^2));
disp(['Iterative ARX-based OE Model RMSE: ', num2str(RMSE_oe)]);

figure;
plot(1:length(y), y, 'b', 1:length(y_oe), y_oe, 'r--');
legend('True Output', 'OE Model Output');
title('OE Model via Iterative ARX Filtering');
xlabel('Time Steps');
ylabel('Output');
grid on;
