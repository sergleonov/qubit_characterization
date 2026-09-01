%% ENGS 53; Extracting T1
% 
% 12/03/2026 Sergei Leonov

close all
clear
clc

data = readtable("time_domain_T1_data.csv");

time = data.time_ns; % nanoseconds
prob = data.p1_meas;

% crop data to extract transient
start = find(prob >= 0.98, 1, 'first');
time = time(start:end); % crop time
time = time - time(1); % start time at zero
prob = prob(start:end); % crop prob

% define the model
exp_model = fittype('A * exp(-x/T1) + C', ...
                   'independent', 'x', 'coefficients', {'A', 'T1', 'C'});

% define initial guess parameters
A = 1;
T1 = 50;
C = 0;

% initial guess
b = [A T1 C];

% extract model parameters
[f, gof, out] = fit(time, prob, exp_model, "StartPoint", b)
A = f.A
T1 = f.T1
C = f.C

% compute predicted model
y = A*exp(-time/T1) + C;

% compute covariance matrix
J = out.Jacobian;
mse = gof.rmse^2;
cov_mat = mse*inv(J'*J);
sqrt(diag(cov_mat))

% plot data and model
fig =figure;
plot(time, prob, '-', 'LineWidth', 1.6);
grid on; hold on;
plot(time, y, '-', 'LineWidth', 2);
legend("Data", "T1 model", "FontSize", 30);
xlabel("Time [ns]", "FontSize", 20);
ylabel("P(|1>)", "FontSize", 20);
set(gca, 'FontSize', 20); 
saveas(fig, "T1.jpeg")
