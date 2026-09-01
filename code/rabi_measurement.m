%% ENGS 53; Extracting Rabi Rate
% 
% 12/03/2026 Sergei Leonov

close all
clear
clc

data = readtable("time_domain_sim.csv");
time = data.time_ns*10^(-3); % convert to microseconds
prob = data.p1_meas;

% define model parameters
A = 1;
omega = 50; % [MHz]
T1 = 0.180; % [us]
C = 0;

% define new model
sin2_fit = fittype('A * sin(4*pi* B * x)^2 * exp(-x/T1)', ...
                   'independent', 'x', 'coefficients', {'A', 'T1', 'B'});
start_pts = [A, T1, omega];

% fit the model
[f, gof, out] = fit(time, prob, sin2_fit, 'StartPoint', start_pts,  "Lower", [0.98 0.180, 0])

% compute model
A = f.A
T1 = f.T1
omega = f.B

y = A*(sin(4*pi*omega*time)).^2.*exp(-time/T1);

% compute covariance matrix
J = out.Jacobian;
mse = gof.rmse^2;
cov_mat = mse*inv(J'*J);
sqrt(diag(cov_mat))

% plot data and model
fig = figure;
plot(time, prob, '-', 'LineWidth', 1.6);
grid on; hold on;
plot(time, y, '-', 'LineWidth', 2);
legend("Data", "Rabi model", "FontSize", 30);
xlabel("Time [us]", "FontSize", 20);
ylabel("P(|1>)", "FontSize", 20);
set(gca, 'FontSize', 20); 
saveas(fig, "Rabi.jpeg")

