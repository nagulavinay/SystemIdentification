 rng(123);
% Sampling time
T0 = 0.1;
t = 0:T0:100;
Sigma = 0.05;
noise = Sigma * randn(size(t));
% Input Signals
step_input = ones(size(t)) + noise;
impulse_input = [1 zeros(1, length(t)-1)] ;
% PRBS generation using prbs.m
n = 10;
m =11;
prbs_input = prbs(n,m) + noise;
prbs_input = prbs_input(1:length(t));
% Noise signals

% Simulate system response for G1
y_step_G1 = lsim(G1, step_input, t) + noise;
y_impulse_G1 = lsim(G1, impulse_input, t) + noise;
y_prbs_G1 = lsim(G1, prbs_input, t) + noise;
% Simulate system response for G2
y_step_G2 = lsim(G2, step_input, t) + noise;
y_impulse_G2 = lsim(G2, impulse_input, t) + noise;
y_impulse_G2_P = lsim(G2, impulse_input, t);
y_prbs_G2 = lsim(G2, prbs_input, t) + noise;
y_prbs_G2_p = lsim(G2, prbs_input, t);
% Response plot for G1
% Plot step response
%figure;
%plot(t, step_input, 'r', t, y_step_G1, 'b');
%legend('Step Input', 'G1 Response');
%title('Step Response (G1)');
% Plot impulse response
%figure;
%plot(t, impulse_input, 'r', t, y_impulse_G1, 'b');
%legend('Impulse Input', 'G1 Response');
%title('Impulse Response (G1)');
% PRBS input signal
%figure;
%plot(t, prbs_input,'k');
%title('PRBS Signal');
%xlabel('Time(s)');
%ylabel('Amplitude');
%grid on; 

% Plot PRBS response
%figure;
%plot(t, prbs_input, 'r', t, y_prbs_G1, 'b');
%legend('PRBS Input', 'G1 Response');
%title('PRBS Response (G1)');

% Response plot for G2
% Plot step response
%figure;
%plot(t, step_input, 'r', t, y_step_G2, 'b');
%legend('Step Input', 'G2 Response');
%title('Step Response (G2)');
% Plot impulse response
%figure;
%plot(t, impulse_input, 'r', t, y_impulse_G2, 'b');
%legend('Impulse Input', 'G2 Response');
%title('Impulse Response (G2)');
% Plot PRBS response
%figure;
%plot(t, prbs_input, 'r', t, y_prbs_G2, 'b');
%legend('PRBS Input', 'G2 Response');
%title('PRBS Response (G2)');