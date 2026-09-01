%% ENGS 53; Extracting Qubit Frequency
% 
% 12/03/2026 Sergei Leonov

close all
clear
clc

data = readtable("frequency_sweep_summary.csv");
freq = data.omega_d_GHz;
prob = data.max_prob_1;

% define the model
lorentz = fittype('A * ((G/2)^2 ./ ((x - wq).^2 + (G/2)^2)) + C', ...
                       'independent', 'x', 'coefficients', {'A', 'G', 'wq', 'C'});

% initial guess
A = 1;
G = 0.5;
wq = 5;
C = min(prob); 
b = [A, G, wq, C];

% fit the model
[f, gof, out] = fit(freq, prob, lorentz, 'StartPoint', b)

% compute the fit
A = f.A
G = f.G
omega_q = f.wq
C = f.C

y = A * ((G/2)^2)./((freq - omega_q).^2 + (G/2)^2) + C;

% covariance matrix
J = out.Jacobian;
mse = gof.rmse^2;
cov_mat = mse*inv(J'*J);
sqrt(diag(cov_mat))

rmse(prob, y)

% plot data and model
fig = figure;
plot(freq, prob, '-', 'LineWidth', 1.6);
grid on; hold on;
plot(freq, y, '-', 'LineWidth', 2);
legend("Data", "Lorentz model", "FontSize", 30);
% title("Frequency Sweep Data Fit")
xlabel("Frequency [GHz]", "FontSize", 20);
ylabel("Maximum P(|1>)", "FontSize", 20);
set(gca, 'FontSize', 20); 
saveas(fig, "freq.jpeg")