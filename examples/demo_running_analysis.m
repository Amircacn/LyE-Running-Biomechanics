%% Demo 4: Real Running Kinematics Analysis
% Full end-to-end processing of 200 Hz optical motion capture running kinematics.
%
% Compares:
%   1. Method 1: Continuous Whole-Sequence Normalization (Mehdizadeh 2014)
%   2. Method 2: Stride-by-Stride Normalization (Warping each stride to 100 pts)
%
% Demonstrates how Method 2 artificially drives long-term LyE into negative values!
%
% Author: LyE Framework Team (2024)

clear; clc; close all;
addpath('../core');

fprintf('====================================================================\n');
fprintf('  REAL RUNNING KINEMATICS LDS ANALYSIS (200 Hz)\n');
fprintf('====================================================================\n\n');

% Load 200 Hz knee kinematics
tbl = readtable('../data/sample_running_knee_angles.csv');
knee_ang = tbl.R_Knee_FlexExt_deg;
fs = 200;
dt = 1 / fs;

% Physical differentiation before normalization
knee_vel = gradient(knee_ang, dt); % deg/s
knee_acc = gradient(knee_vel, dt); % deg/s^2

% Detect strides
[strides, cinfo] = detect_gait_cycles(knee_ang, fs, 'MinStrideTime', 0.5);
fprintf('Detected %d strides (Mean stride time: %.3f s, Cadence: %.1f spm)\n', ...
        cinfo.n_strides, cinfo.mean_stride_time, cinfo.cadence_spm);

n_use = min(80, cinfo.n_strides);
pts_total = n_use * 100;
strides_use = strides(1:(n_use+1));

% -------------------------------------------------------------------------
% Normalization Method 1: Continuous Whole-Sequence
% -------------------------------------------------------------------------
sig_crop = knee_vel(strides_use(1):strides_use(end));
t_crop = 1:length(sig_crop);
vel_m1 = interp1(t_crop, sig_crop, linspace(1, t_crop(end), pts_total), 'spline')';

% -------------------------------------------------------------------------
% Normalization Method 2: Stride-by-Stride Warping
% -------------------------------------------------------------------------
vel_m2 = zeros(pts_total, 1);
for s = 1:n_use
    st_pts = sig_crop(strides_use(s)-strides_use(1)+1 : strides_use(s+1)-strides_use(1)+1);
    t_st = 1:length(st_pts);
    vel_m2(((s-1)*100+1):(s*100)) = interp1(t_st, st_pts, linspace(1, t_st(end), 100), 'spline')';
end

% -------------------------------------------------------------------------
% Estimate LyE for Both Methods
% -------------------------------------------------------------------------
fprintf('\nCalculating LDS for Method 1 (Continuous Whole-Sequence)...\n');
res_m1 = estimate_rosenstein_lye(vel_m1, ...
    'EmbeddingDimension', 5, 'EmbeddingDelay', 10, 'TheilerWindow', 50, ...
    'NumNeighbors', 15, 'DivergenceWindow', 1000, 'SamplesPerStride', 100, ...
    'ShortTermRange', [0 0.5], 'LongTermRange', [4 10]);

fprintf('Calculating LDS for Method 2 (Stride-by-Stride Warping)...\n');
res_m2 = estimate_rosenstein_lye(vel_m2, ...
    'EmbeddingDimension', 5, 'EmbeddingDelay', 10, 'TheilerWindow', 50, ...
    'NumNeighbors', 15, 'DivergenceWindow', 1000, 'SamplesPerStride', 100, ...
    'ShortTermRange', [0 0.5], 'LongTermRange', [4 10]);

% -------------------------------------------------------------------------
% Results & Plotting
% -------------------------------------------------------------------------
fprintf('\n====================================================================\n');
fprintf('Condition                       | LyE_short (per stride) | LyE_long (per stride) | R2 (short)\n');
fprintf('--------------------------------|------------------------|-----------------------|-----------\n');
fprintf('Method 1 (Whole-Sequence)       | %22.4f | %21.4f | %10.3f\n', ...
    res_m1.LyE_short, res_m1.LyE_long, res_m1.fit_short.R2);
fprintf('Method 2 (Stride-by-Stride)     | %22.4f | %21.4f | %10.3f\n', ...
    res_m2.LyE_short, res_m2.LyE_long, res_m2.fit_short.R2);
fprintf('====================================================================\n\n');

figure('Name', 'Real Running LDS Comparison', 'Position', [100 100 950 450], 'Color', 'w');
t_str = (1:1000) / 100;

subplot(1, 2, 1);
plot(t_str, res_m1.div_curve, 'b-', 'LineWidth', 1.8); hold on;
plot(t_str, res_m2.div_curve, 'r--', 'LineWidth', 1.8);
xlabel('Time (strides)'); ylabel('\langle ln(divergence) \rangle');
title('Running Knee Velocity Divergence Curves');
legend('Method 1: Whole-Sequence (Continuous)', 'Method 2: Stride-by-Stride (Warped)', 'Location', 'southeast');
grid on;

subplot(1, 2, 2);
bar_data = [res_m1.LyE_short, res_m2.LyE_short; res_m1.LyE_long, res_m2.LyE_long];
bar(bar_data);
set(gca, 'XTickLabel', {'\lambda_S (0 - 0.5 strides)', '\lambda_L (4 - 10 strides)'});
ylabel('Lyapunov Exponent (per stride)');
title('Impact of Normalization on LyE Exponents');
legend('Method 1', 'Method 2', 'Location', 'northwest');
grid on;
