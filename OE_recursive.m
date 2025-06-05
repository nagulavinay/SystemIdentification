% === OE MODEL from Recursive Filtering ===
% Step 1: Get F(q) from ARX estimate (same F as before)
F_oe = [theta(1:na).'];  % include leading 1
if any(abs(roots(F_oe)) >= 1)
    error('Unstable filter F(q) detected!');
end

% Step 2: Apply 1/F(q) filter to both y and u
u = prbs_input;
y = y_prbs_G2;
y_oe_filt = filter(1, F_oe, y);
u_oe_filt = filter(1, F_oe, u);

% Step 3: Create regression matrix (no AR part in OE)
X_oe = zeros(length(y_oe_filt), nb);
for i = 1:length(y_oe_filt)
    for j = 1:nb
        if (i - nk - j + 1) > 0
            X_oe(i, j) = u_oe_filt(i - nk - j + 1);
        end
    end
end

% Step 4: Estimate OE model parameters (only B coefficients)
theta_oe = (X_oe' * X_oe) \ (X_oe' * y_oe_filt);

% Step 5: Simulate OE output using estimated F and B
y_oe = zeros(size(y));
for i = 1:length(y)
    for j = 1:nb
        if (i - nk - j + 1) > 0
            y_oe(i) = y_oe(i) + theta_oe(j) * u(i - nk - j + 1);
        end
    end
    for j = 1:na
        if i - j > 0
            y_oe(i) = y_oe(i) - F_oe(j) * y_oe(i - j);
        end
    end
end

% Step 6: Compute RMSE between OE output and true output
N_oe = min(length(y_oe_filt), length(y_oe));
y_filt_trim = y_oe_filt(1:N_min);
y_oe_filt_trim = y_oe(1:N_min);
RMSE_oe = sqrt(mean((y_filt_trim(:) - y_oe_filt_trim(:)).^2));
disp(['Recursive Filtering OE Model RMSE: ', num2str(RMSE_oe)]);

% Plot OE vs. true
figure;
plot(1:N_oe, y(1:N_oe), 'b', 'LineWidth', 1.5); hold on;
plot(1:N_oe, y_oe(1:N_oe), 'g--', 'LineWidth', 1.5); hold off;
legend('True Output', 'Estimated OE Model');
title('OE Model Output vs. True Output');
xlabel('Time Steps');
ylabel('Output');
grid on;
