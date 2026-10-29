%% Demo 1: Quickstart Minimal Example
% Minimal, clear example of calculating Local Dynamic Stability (LDS / LyE)
% from running knee joint kinematics using the unified Rosenstein framework.
%
% Author: LyE Framework Team (2024)

clear; clc; close all;
addpath('../core');

% 1. Load sample 200 Hz running knee kinematics
data_file = '../data/sample_running_knee_angles.csv';
tbl = readtable(data_file);
knee_ang = tbl.R_Knee_FlexExt_deg;
fs_raw = 200; % 200 Hz

% 2. Differentiate at raw 200 Hz to get angular velocity (deg/s)
knee_vel = gradient(knee_ang, 1/fs_raw);

% 3. Automated stride detection
[stride_locs, info] = detect_gait_cycles(knee_ang, fs_raw, 'MinStrideTime', 0.5);
fprintf('Detected %d strides (Cadence: %.1f spm)\n', info.n_strides, info.cadence_spm);

% 4. Continuous whole-sequence normalization (80 strides -> 8000 points)
n_strides = min(80, info.n_strides);
pts_total = n_strides * 100;
sig_crop = knee_vel(stride_locs(1):stride_locs(n_strides+1));
norm_vel = interp1(1:length(sig_crop), sig_crop, linspace(1, length(sig_crop), pts_total), 'spline')';

% 5. Estimate Largest Lyapunov Exponent (LDS)
results = estimate_rosenstein_lye(norm_vel, ...
    'EmbeddingDimension', 5, ...
    'EmbeddingDelay', 10, ...
    'TheilerWindow', 50, ...
    'NumNeighbors', 15, ...
    'DivergenceWindow', 1000, ...
    'SamplesPerStride', 100, ...
    'ShortTermRange', [0 0.5], ...
    'LongTermRange', [4 10], ...
    'TimeUnit', 'strides', ...
    'PlotResults', true);

% 6. Display results
fprintf('\n--- RUNNING KINEMATICS LDS RESULTS ---\n');
fprintf('Short-term LyE (lambda_S): %.4f %s (R^2 = %.3f)\n', ...
    results.LyE_short, results.LyE_short_unit, results.fit_short.R2);
fprintf('Long-term LyE  (lambda_L): %.4f %s\n', ...
    results.LyE_long, results.LyE_long_unit);
fprintf('Execution Time          : %.3f seconds\n', results.execution_time);
