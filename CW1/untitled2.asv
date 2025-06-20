% (iii) Additive perturbation analysis with improved numerical evaluation

% Define frequency range
omega = logspace(-1, 2, 500); % Frequency range from 0.1 to 100 rad/s

cp_range = [0.925000000000000	1.07500000000000];
cv_range = [0.900000000000000	1.10000000000000];

% Compute perturbed G(s) for extreme uncertainty cases
G_low = 1 / (m * s^2 + cv_range(1) * s + cp_range(1)); % Lower bound system
G_high = 1 / (m * s^2 + cv_range(2) * s + cp_range(2)); % Upper bound system

% Compute frequency responses
[mag_nom, ~] = bode(G_nom, omega);
[mag_low, ~] = bode(G_low, omega);
[mag_high, ~] = bode(G_high, omega);

% Convert cell arrays from bode() output to matrices
mag_nom = squeeze(mag_nom);
mag_low = squeeze(mag_low);
mag_high = squeeze(mag_high);

% Compute the worst-case perturbation bound: b(ω) = max(|G_high(jω) - G_nom(jω)|, |G_nom(jω) - G_low(jω)|)
Delta_mag = max(abs(mag_high - mag_nom), abs(mag_nom - mag_low));

% Bode magnitude plot for perturbation bound b(ω)
figure;
semilogx(omega, 20*log10(Delta_mag), 'r', 'LineWidth', 1.5);
grid on;
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');
title('Bound on Additive Perturbation b(ω)');
legend('|Δ(jω)|');