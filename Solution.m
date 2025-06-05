s = tf('s');
T0 = 0.1;
t = 0:T0:100;
% System Definition
G1 = 10 / (s^3 + 5.8*s^2 + 8*s + 20);
G2 = (0.25 / (s^2 + 1.5*s + 0.5)) * exp(-s);
% Step Response
%figure; step(G1); title('Step Response of G1')
%figure; step(G2); title ('Step Response of G2');
% Impulse Response
%figure; impulse(G1); title('Impulse Response of G1');
%figure; impulse(G2); title ('Impulse Response of G2');
% Bode Plot
%figure; bode(G1); title('Bode Plot of G1');
%figure; bode(G2); title('Bode Plot of G2');
% Pole Locations
poles_G1 = pole(G1);
disp('Poles of G1');
disp(poles_G1);
poles_G2 = pole(G2);
disp('Poles of G2');
disp(poles_G2);