function compare_algorithms_and_tune_user_data()
%COMPARE_ALGORITHMS_AND_TUNE_USER_DATA Head-to-head empirical comparison of
% raw Mehdizadeh code, raw NONAN (LyE_R) code, and our new framework directly
% on the user's running knee kinematics dataset, followed by systematic tuning.
%
% Experiments:
%   1. Direct Head-to-Head Comparison on user's running data:
%      - Original lds_calc_mehdizadeh.m (k=1 and k=15)
%      - Original LyE_R.m (NONAN: k=1, Theiler = 0.8*tau)
%      - Our Unified Framework (estimate_rosenstein_lye)
%   2. Systematic sweep of k (1 to 30) on user's joint angle & velocity.
%   3. Systematic sweep of Theiler Window:
%      - NONAN rule: W = 0.8*tau = 8 samples
%      - Mehdizadeh rule: W = 0.5*period = 50 samples
%      - Full-stride rule: W = 1.0*period = 100 samples
%      - Fine grid: W in [4, 8, 10, 20, 30, 50, 75, 100, 150]
%   4. Embedding delay (tau) and dimension (d_E) sweeps.
%
% Author: LyE Framework Team (2024)

    fprintf('====================================================================\n');
    fprintf('  HEAD-TO-HEAD COMPARISON & PARAMETER TUNING ON USER RUNNING DATA\n');
    fprintf('====================================================================\n\n');

    base_dir = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(base_dir, 'core'));
    addpath(fullfile(base_dir, 'original_codes'));

    % Load user's running data: Subject 01, Control condition
    data_file = 'c:\Users\Amir\Downloads\LyPE\Data\sub01\kinematics\control\R_knee_angles.csv';
    if ~exist(data_file, 'file')
        error('Data file not found: %s', data_file);
    end

    tbl = readtable(data_file);
    raw_angle = tbl.R_Knee_FlexExt_deg;
    fs_raw = 200; % 200 Hz
    dt_raw = 1 / fs_raw;

    % Physical differentiation on raw time axis
    raw_vel = gradient(raw_angle, dt_raw); % deg/s

    % Detect gait strides
    [locs, ~] = detect_gait_cycles(raw_angle, fs_raw, 'MinStrideTime', 0.5);
    n_strides = min(80, length(locs) - 1);
    target_pts = n_strides * 100;
    strides_use = locs(1:(n_strides + 1));

    % Continuous Whole-Sequence Normalization (Method 1)
    sig_crop_ang = raw_angle(strides_use(1):strides_use(end));
    sig_crop_vel = raw_vel(strides_use(1):strides_use(end));
    t_orig = 1:length(sig_crop_ang);
    t_norm = linspace(1, t_orig(end), target_pts);

    ang_norm = interp1(t_orig, sig_crop_ang, t_norm, 'spline')';
    vel_norm = interp1(t_orig, sig_crop_vel, t_norm, 'spline')';

    fprintf('User Data Loaded: 80 strides, %d resampled points (100 samples/stride)\n\n', target_pts);

    % =====================================================================
    % EXPERIMENT 1: DIRECT HEAD-TO-HEAD ON UNTOUCHED ORIGINAL CODES
    % =====================================================================
    fprintf('--------------------------------------------------------------------\n');
    fprintf('EXPERIMENT 1: Running Untouched Original Codes vs Our Framework\n');
    fprintf('--------------------------------------------------------------------\n');

    dim = 5;
    tau = 10;
    period = 1; % 1 stride in normalized time
    fs_norm = 100; % 100 samples per stride
    ws = 10; % 10 strides window

    % --- Code 1A: Raw lds_calc_mehdizadeh (k = 1) ---
    [S_meh, ~] = embed_state_space(vel_norm, dim, tau);
    [div_meh_k1, lds_meh_k1] = lds_calc_mehdizadeh(S_meh, ws, fs_norm, period, 1, 0);
    % Note: lds_calc_mehdizadeh fits against seconds: 1/fs:1/fs:...
    % To convert to per stride, multiply by stride time (1 second here since fs_norm=100 & period=1)
    lye_meh_k1_s = lds_meh_k1(1);
    lye_meh_k1_l = lds_meh_k1(2);

    % --- Code 1B: Raw lds_calc_mehdizadeh (k = 15) ---
    [div_meh_k15, lds_meh_k15] = lds_calc_mehdizadeh(S_meh, ws, fs_norm, period, 15, 0);
    lye_meh_k15_s = lds_meh_k15(1);
    lye_meh_k15_l = lds_meh_k15(2);

    % --- Code 2: Raw LyE_R.m (NONAN toolbox) ---
    % LyE_R(X, Fs, tau, dim, slope, MeanPeriod, file)
    % LyE_R defines exclusion window as: range_exclude = i - round(tau*0.8) : round(i + tau*0.8)
    % slope = [0.05 0.5 4 10] strides, MeanPeriod = 100 samples, Fs = 100 Hz
    [lyes_nonan, lyel_nonan, out_nonan] = LyE_R(vel_norm, fs_norm, tau, dim, [0.05 0.5 4 10], 1.0, 0);

    % --- Code 3A: Our Unified Framework (k = 1, Theiler = 50 samples = 0.5 stride) ---
    res_our_k1 = estimate_rosenstein_lye(vel_norm, ...
        'EmbeddingDimension', dim, 'EmbeddingDelay', tau, 'TheilerWindow', 50, ...
        'NumNeighbors', 1, 'DivergenceWindow', 1000, 'SamplesPerStride', 100, ...
        'ShortTermRange', [0 0.5], 'LongTermRange', [4 10]);

    % --- Code 3B: Our Unified Framework (k = 15, Theiler = 50 samples = 0.5 stride) ---
    res_our_k15 = estimate_rosenstein_lye(vel_norm, ...
        'EmbeddingDimension', dim, 'EmbeddingDelay', tau, 'TheilerWindow', 50, ...
        'NumNeighbors', 15, 'DivergenceWindow', 1000, 'SamplesPerStride', 100, ...
        'ShortTermRange', [0 0.5], 'LongTermRange', [4 10]);

    % --- Code 3C: Our Framework with NONAN Theiler Window (Theiler = 8 samples = 0.8*tau) ---
    res_our_nonan_theiler = estimate_rosenstein_lye(vel_norm, ...
        'EmbeddingDimension', dim, 'EmbeddingDelay', tau, 'TheilerWindow', round(0.8*tau), ...
        'NumNeighbors', 1, 'DivergenceWindow', 1000, 'SamplesPerStride', 100, ...
        'ShortTermRange', [0 0.5], 'LongTermRange', [4 10]);

    fprintf('\n>>> EXPERIMENT 1 RESULTS (Knee Angular Velocity):\n');
    fprintf('1. Raw Mehdizadeh (k=1)             : LyE_S = %.4f | LyE_L = %.4f\n', lye_meh_k1_s, lye_meh_k1_l);
    fprintf('2. Raw Mehdizadeh (k=15)            : LyE_S = %.4f | LyE_L = %.4f\n', lye_meh_k15_s, lye_meh_k15_l);
    fprintf('3. Raw NONAN LyE_R (k=1, W=0.8*tau) : LyE_S = %.4f | LyE_L = %.4f\n', lyes_nonan, lyel_nonan);
    fprintf('4. Our Framework (k=1, W=50 samples): LyE_S = %.4f | LyE_L = %.4f | R2 = %.3f\n', ...
            res_our_k1.LyE_short, res_our_k1.LyE_long, res_our_k1.fit_short.R2);
    fprintf('5. Our Framework (k=15, W=50 sample): LyE_S = %.4f | LyE_L = %.4f | R2 = %.3f\n', ...
            res_our_k15.LyE_short, res_our_k15.LyE_long, res_our_k15.fit_short.R2);
    fprintf('6. Our Framework with NONAN Theiler : LyE_S = %.4f | LyE_L = %.4f | R2 = %.3f\n\n', ...
            res_our_nonan_theiler.LyE_short, res_our_nonan_theiler.LyE_long, res_our_nonan_theiler.fit_short.R2);

    % =====================================================================
    % EXPERIMENT 2: SYSTEMATIC SWEEP OF NEIGHBORS (k = 1 to 30)
    % =====================================================================
    fprintf('--------------------------------------------------------------------\n');
    fprintf('EXPERIMENT 2: Sweeping k (1 to 30) on User Knee Velocity & Angle\n');
    fprintf('--------------------------------------------------------------------\n');

    k_grid = [1, 2, 3, 5, 8, 10, 12, 15, 18, 20, 25, 30];
    lye_k_vel = zeros(size(k_grid));
    lye_k_ang = zeros(size(k_grid));
    r2_k_vel  = zeros(size(k_grid));
    r2_k_ang  = zeros(size(k_grid));

    for i = 1:length(k_grid)
        k_val = k_grid(i);

        % On Knee Velocity
        r_v = estimate_rosenstein_lye(vel_norm, ...
            'EmbeddingDimension', 5, 'EmbeddingDelay', 10, 'TheilerWindow', 50, ...
            'NumNeighbors', k_val, 'DivergenceWindow', 1000, 'SamplesPerStride', 100, ...
            'ShortTermRange', [0 0.5], 'LongTermRange', [4 10]);
        lye_k_vel(i) = r_v.LyE_short;
        r2_k_vel(i)  = r_v.fit_short.R2;

        % On Knee Angle
        r_a = estimate_rosenstein_lye(ang_norm, ...
            'EmbeddingDimension', 5, 'EmbeddingDelay', 10, 'TheilerWindow', 50, ...
            'NumNeighbors', k_val, 'DivergenceWindow', 1000, 'SamplesPerStride', 100, ...
            'ShortTermRange', [0 0.5], 'LongTermRange', [4 10]);
        lye_k_ang(i) = r_a.LyE_short;
        r2_k_ang(i)  = r_a.fit_short.R2;

        fprintf('  k = %2d | Knee Velocity LyE_S = %.4f (R2=%.3f) | Knee Angle LyE_S = %.4f (R2=%.3f)\n', ...
                k_val, lye_k_vel(i), r2_k_vel(i), lye_k_ang(i), r2_k_ang(i));
    end

    % Check derivative of LyE w.r.t k to detect where plateau forms (|dLyE/dk| < threshold)
    dk_vel = abs(diff(lye_k_vel) ./ diff(k_grid));
    fprintf('\nPlateau check on Knee Velocity: d(LyE)/dk between successive k steps:\n');
    for i = 1:length(dk_vel)
        fprintf('  Between k=%2d and k=%2d : rate of change = %.5f\n', k_grid(i), k_grid(i+1), dk_vel(i));
    end

    % =====================================================================
    % EXPERIMENT 3: SYSTEMATIC SWEEP OF THEILER WINDOW (W)
    % =====================================================================
    fprintf('\n--------------------------------------------------------------------\n');
    fprintf('EXPERIMENT 3: Theiler Window Comparison (NONAN vs Mehdizadeh vs Fine Grid)\n');
    fprintf('--------------------------------------------------------------------\n');

    theiler_grid = [4, 8, 10, 15, 20, 30, 50, 75, 100, 150];
    lye_th_vel = zeros(size(theiler_grid));
    lye_th_ang = zeros(size(theiler_grid));

    for i = 1:length(theiler_grid)
        th_val = theiler_grid(i);

        r_v = estimate_rosenstein_lye(vel_norm, ...
            'EmbeddingDimension', 5, 'EmbeddingDelay', 10, 'TheilerWindow', th_val, ...
            'NumNeighbors', 15, 'DivergenceWindow', 1000, 'SamplesPerStride', 100, ...
            'ShortTermRange', [0 0.5], 'LongTermRange', [4 10]);
        lye_th_vel(i) = r_v.LyE_short;

        r_a = estimate_rosenstein_lye(ang_norm, ...
            'EmbeddingDimension', 5, 'EmbeddingDelay', 10, 'TheilerWindow', th_val, ...
            'NumNeighbors', 15, 'DivergenceWindow', 1000, 'SamplesPerStride', 100, ...
            'ShortTermRange', [0 0.5], 'LongTermRange', [4 10]);
        lye_th_ang(i) = r_a.LyE_short;

        fprintf('  Theiler W = %3d samples (%4.1f%% stride) | Velocity LyE_S = %.4f | Angle LyE_S = %.4f\n', ...
                th_val, th_val/100*100, lye_th_vel(i), lye_th_ang(i));
    end

    % =====================================================================
    % EXPERIMENT 4: EMBEDDING DELAY (tau) AND DIMENSION (d_E) SWEEPS
    % =====================================================================
    fprintf('\n--------------------------------------------------------------------\n');
    fprintf('EXPERIMENT 4: Embedding Delay and Dimension Sweeps\n');
    fprintf('--------------------------------------------------------------------\n');

    delays = [5, 8, 10, 12, 15, 20];
    lye_delays = zeros(size(delays));
    for i = 1:length(delays)
        r = estimate_rosenstein_lye(vel_norm, ...
            'EmbeddingDimension', 5, 'EmbeddingDelay', delays(i), 'TheilerWindow', 50, ...
            'NumNeighbors', 15, 'DivergenceWindow', 1000, 'SamplesPerStride', 100, ...
            'ShortTermRange', [0 0.5], 'LongTermRange', [4 10]);
        lye_delays(i) = r.LyE_short;
        fprintf('  Delay tau = %2d samples | Knee Velocity LyE_S = %.4f (R2=%.3f)\n', ...
                delays(i), lye_delays(i), r.fit_short.R2);
    end

    dims = [3, 4, 5, 6, 7];
    lye_dims = zeros(size(dims));
    for i = 1:length(dims)
        r = estimate_rosenstein_lye(vel_norm, ...
            'EmbeddingDimension', dims(i), 'EmbeddingDelay', 10, 'TheilerWindow', 50, ...
            'NumNeighbors', 15, 'DivergenceWindow', 1000, 'SamplesPerStride', 100, ...
            'ShortTermRange', [0 0.5], 'LongTermRange', [4 10]);
        lye_dims(i) = r.LyE_short;
        fprintf('  Dimension d_E = %d | Knee Velocity LyE_S = %.4f (R2=%.3f)\n', ...
                dims(i), lye_dims(i), r.fit_short.R2);
    end

    % =====================================================================
    % SAVE RESULTS & PLOT FIGURES
    % =====================================================================
    results_struct.exp1.meh_k1_s = lye_meh_k1_s;
    results_struct.exp1.meh_k1_l = lye_meh_k1_l;
    results_struct.exp1.meh_k15_s = lye_meh_k15_s;
    results_struct.exp1.meh_k15_l = lye_meh_k15_l;
    results_struct.exp1.nonan_s = lyes_nonan;
    results_struct.exp1.nonan_l = lyel_nonan;
    results_struct.exp1.our_k1 = res_our_k1;
    results_struct.exp1.our_k15 = res_our_k15;
    results_struct.exp1.our_nonan_theiler = res_our_nonan_theiler;

    results_struct.exp2.k_grid = k_grid;
    results_struct.exp2.lye_k_vel = lye_k_vel;
    results_struct.exp2.lye_k_ang = lye_k_ang;
    results_struct.exp2.r2_k_vel = r2_k_vel;

    results_struct.exp3.theiler_grid = theiler_grid;
    results_struct.exp3.lye_th_vel = lye_th_vel;
    results_struct.exp3.lye_th_ang = lye_th_ang;

    results_struct.exp4.delays = delays;
    results_struct.exp4.lye_delays = lye_delays;
    results_struct.exp4.dims = dims;
    results_struct.exp4.lye_dims = lye_dims;

    res_dir = fullfile(base_dir, 'results');
    save(fullfile(res_dir, 'user_data_comparison_and_tuning.mat'), 'results_struct', '-v7');

    % Generate Dedicated Figure for Tuning
    fig = figure('Name', 'UserData_Comparison_and_Tuning', 'Position', [100 100 1100 700], 'Color', 'w', 'Visible', 'off');

    % Subplot 1: Direct Method Comparison
    subplot(2, 2, 1);
    bar_comp = [lye_meh_k1_s, lyes_nonan, res_our_k15.LyE_short];
    b = bar(bar_comp, 'FaceColor', [0.3 0.6 0.8]);
    set(gca, 'XTickLabel', {'Mehdizadeh (k=1)', 'NONAN (W=0.8\tau)', 'Our Framework (k=15)'});
    ylabel('Short-term LyE \lambda_S (per stride)');
    title('Comparison of Algorithms on User Knee Velocity');
    grid on;

    % Subplot 2: k sweep (Plateau demonstration)
    subplot(2, 2, 2);
    plot(k_grid, lye_k_vel, 'b-o', 'LineWidth', 2); hold on;
    plot(k_grid, lye_k_ang, 'r-s', 'LineWidth', 1.8);
    xline(15, 'g--', 'Optimal k = 15 Plateau', 'LineWidth', 1.5);
    xlabel('Number of Nearest Neighbors (k)');
    ylabel('Estimated LyE \lambda_S (per stride)');
    title('Neighbor Sweep: Plateau Formation');
    legend('Knee Angular Velocity', 'Knee Angle', 'k=15 Plateau', 'Location', 'best');
    grid on;

    % Subplot 3: Theiler window sweep (NONAN vs Mehdizadeh vs Optimal)
    subplot(2, 2, 3);
    plot(theiler_grid, lye_th_vel, 'b-d', 'LineWidth', 2); hold on;
    xline(8, 'r:', 'NONAN Rule: W = 0.8\tau (8 samples)', 'LineWidth', 1.5);
    xline(50, 'g--', 'Mehdizadeh Rule: W = 0.5 Stride (50 samples)', 'LineWidth', 1.5);
    xlabel('Theiler Window W (samples)');
    ylabel('Estimated LyE \lambda_S (per stride)');
    title('Theiler Window Sensitivity on Running Data');
    legend('Knee Velocity LyE', 'NONAN W=8', '0.5 Stride W=50', 'Location', 'best');
    grid on;

    % Subplot 4: Delay sensitivity
    subplot(2, 2, 4);
    plot(delays, lye_delays, 'm-^', 'LineWidth', 2);
    xlabel('Embedding Delay \tau (samples)');
    ylabel('Estimated LyE \lambda_S (per stride)');
    title('Embedding Delay Sensitivity (\tau)');
    xline(10, 'k--', 'Recommended \tau = 10');
    grid on;

    fig_path = fullfile(res_dir, 'figures', 'UserData_Comparison_and_Tuning.png');
    saveas(fig, fig_path);
    close(fig);
    fprintf('\nTuning figure saved to: %s\n', fig_path);
    fprintf('Results saved to: %s\n\n', fullfile(res_dir, 'user_data_comparison_and_tuning.mat'));

end
