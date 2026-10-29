%% Demo 2: Center of Mass (CoM) Velocity vs Joint Angle Dynamics
% Demonstrates why Center of Mass velocity and Joint Angle kinematics
% exhibit distinct phase space topologies and require tailored parameter tuning.
%
% Key Takeaway:
%   - CoM Velocity reflects whole-body linear momentum balance (LyE ~ 0.7 - 0.9).
%   - Joint Flexion/Extension Angles reflect segment articulation dynamics (LyE ~ 1.5 - 2.0).
%   - Neither number is "wrong"; each quantifies dynamic sensitivity at a different
%     hierarchical level of the neuromusculoskeletal system!
%
% Author: LyE Framework Team (2024)

clear; clc; close all;
addpath('../core');

fprintf('====================================================================\n');
fprintf('  DEMO: CoM VELOCITY VS JOINT ANGLE KINEMATICS DYNAMICS\n');
fprintf('====================================================================\n\n');

% -------------------------------------------------------------------------
% 1. Load CoM Velocity Data
% -------------------------------------------------------------------------
com_file = '../data/sample_com_velocity.mat';
if ~exist(com_file, 'file')
    error('File not found: %s. Run from examples/ folder.', com_file);
end
loaded_com = load(com_file);
com_vel = loaded_com.CoM_vel;
lhs = loaded_com.lhs_events;
n_strides_com = min(80, length(lhs) - 1);

% Whole-sequence continuous normalization to 100 samples/stride
pts_com = n_strides_com * 100;
sig_crop_com = com_vel(lhs(1):lhs(n_strides_com+1));
t_orig_c = 1:length(sig_crop_com);
com_vel_norm = interp1(t_orig_c, sig_crop_com, linspace(1, t_orig_c(end), pts_com), 'spline')';

% -------------------------------------------------------------------------
% 2. Load Knee Joint Angle Data (200 Hz raw)
% -------------------------------------------------------------------------
knee_file = '../data/sample_running_knee_angles.csv';
tbl_knee = readtable(knee_file);
knee_ang = tbl_knee.R_Knee_FlexExt_deg;
fs_knee = 200;
dt_knee = 1 / fs_knee;

% Compute physical angular velocity at raw 200 Hz
knee_vel = gradient(knee_ang, dt_knee);

% Detect strides on knee flexion
[knee_strides, ~] = detect_gait_cycles(knee_ang, fs_knee, 'MinStrideTime', 0.5);
n_strides_knee = min(80, length(knee_strides) - 1);
pts_knee = n_strides_knee * 100;

sig_crop_k_ang = knee_ang(knee_strides(1):knee_strides(n_strides_knee+1));
sig_crop_k_vel = knee_vel(knee_strides(1):knee_strides(n_strides_knee+1));
t_orig_k = 1:length(sig_crop_k_ang);

knee_ang_norm = interp1(t_orig_k, sig_crop_k_ang, linspace(1, t_orig_k(end), pts_knee), 'spline')';
knee_vel_norm = interp1(t_orig_k, sig_crop_k_vel, linspace(1, t_orig_k(end), pts_knee), 'spline')';

% -------------------------------------------------------------------------
% 3. Calculate LyE for All Three Biomechanical Variables
% -------------------------------------------------------------------------
fprintf('Estimating LDS for CoM Medio-Lateral Velocity (k=15)...\n');
res_com = estimate_rosenstein_lye(com_vel_norm, ...
    'EmbeddingDimension', 5, 'EmbeddingDelay', 10, 'TheilerWindow', 50, ...
    'NumNeighbors', 15, 'DivergenceWindow', 1000, 'SamplesPerStride', 100, ...
    'ShortTermRange', [0 0.5], 'LongTermRange', [4 10], 'TimeUnit', 'strides');

fprintf('Estimating LDS for Knee Joint Angle (k=15)...\n');
res_knee_ang = estimate_rosenstein_lye(knee_ang_norm, ...
    'EmbeddingDimension', 5, 'EmbeddingDelay', 10, 'TheilerWindow', 50, ...
    'NumNeighbors', 15, 'DivergenceWindow', 1000, 'SamplesPerStride', 100, ...
    'ShortTermRange', [0 0.5], 'LongTermRange', [4 10], 'TimeUnit', 'strides');

fprintf('Estimating LDS for Knee Joint Angular Velocity (k=15)...\n');
res_knee_vel = estimate_rosenstein_lye(knee_vel_norm, ...
    'EmbeddingDimension', 5, 'EmbeddingDelay', 10, 'TheilerWindow', 50, ...
    'NumNeighbors', 15, 'DivergenceWindow', 1000, 'SamplesPerStride', 100, ...
    'ShortTermRange', [0 0.5], 'LongTermRange', [4 10], 'TimeUnit', 'strides');

% -------------------------------------------------------------------------
% 4. Comparative Results Table
% -------------------------------------------------------------------------
fprintf('\n====================================================================\n');
fprintf('  BIOMECHANICAL COMPARISON: CoM vs JOINT KINEMATICS\n');
fprintf('====================================================================\n');
fprintf('Variable                       | LyE_short (per stride) | LyE_long (per stride) | R2 (short)\n');
fprintf('-------------------------------|------------------------|-----------------------|-----------\n');
fprintf('CoM Medio-Lateral Velocity     | %22.4f | %21.4f | %10.3f\n', ...
    res_com.LyE_short, res_com.LyE_long, res_com.fit_short.R2);
fprintf('Knee Flexion/Extension Angle   | %22.4f | %21.4f | %10.3f\n', ...
    res_knee_ang.LyE_short, res_knee_ang.LyE_long, res_knee_ang.fit_short.R2);
fprintf('Knee Angular Velocity (dTheta) | %22.4f | %21.4f | %10.3f\n', ...
    res_knee_vel.LyE_short, res_knee_vel.LyE_long, res_knee_vel.fit_short.R2);
fprintf('====================================================================\n\n');

% -------------------------------------------------------------------------
% 5. Visualization
% -------------------------------------------------------------------------
figure('Name', 'CoM vs Joint Angle Dynamics', 'Position', [100 100 1000 450], 'Color', 'w');
t_str = (1:1000) / 100;

subplot(1, 2, 1);
plot(t_str, res_com.div_curve(1:1000), 'b-', 'LineWidth', 1.8); hold on;
plot(t_str, res_knee_vel.div_curve(1:1000), 'r--', 'LineWidth', 1.8);
plot(t_str, res_knee_ang.div_curve(1:1000), 'k-.', 'LineWidth', 1.5);
xlabel('Time (strides)'); ylabel('\langle ln(divergence) \rangle');
title('Logarithmic Divergence Curves');
legend('CoM ML Velocity', 'Knee Angular Velocity', 'Knee Flexion Angle', 'Location', 'southeast');
grid on;

subplot(1, 2, 2);
bar_vals = [res_com.LyE_short, res_knee_vel.LyE_short, res_knee_ang.LyE_short];
bar(bar_vals, 'FaceColor', [0.2 0.5 0.8]);
set(gca, 'XTickLabel', {'CoM Velocity', 'Knee Ang. Vel.', 'Knee Angle'});
ylabel('Short-term LyE \lambda_S (per stride)');
title('Comparison of Local Dynamic Stability Exponents');
grid on;
