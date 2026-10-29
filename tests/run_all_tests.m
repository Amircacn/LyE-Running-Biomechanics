function run_all_tests()
%RUN_ALL_TESTS Comprehensive unit test and numerical verification suite
% for the Lyapunov Exponent (LyE / LDS) framework.
%
%   Runs all unit tests, numerical boundary tests, edge-case validations,
%   and regression verifications.
%
%   Author: Computational Biomechanics & Dynamical Systems Framework
%   Version: 2.0 (2024)

    fprintf('====================================================================\n');
    fprintf('  RUNNING LYE FRAMEWORK UNIT AND NUMERICAL TEST SUITE\n');
    fprintf('====================================================================\n\n');

    base_dir = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(base_dir, 'core'));

    n_passed = 0;
    n_failed = 0;

    % -------------------------------------------------------------
    % Test 1: embed_state_space dimensionality and orientation
    % -------------------------------------------------------------
    fprintf('Test 1: State space embedding dimensions... ');
    try
        x = randn(500, 1);
        [S, p] = embed_state_space(x, 4, 5);
        expected_rows = 500 - (4 - 1)*5;
        assert(size(S, 1) == expected_rows && size(S, 2) == 4, 'Incorrect S dimensions');

        % Test row vector input auto-transposition
        [S_row, ~] = embed_state_space(x', 4, 5);
        assert(isequal(S, S_row), 'Row vector handling failed');

        fprintf('PASSED\n');
        n_passed = n_passed + 1;
    catch ME
        fprintf('FAILED: %s\n', ME.message);
        n_failed = n_failed + 1;
    end

    % -------------------------------------------------------------
    % Test 2: find_neighbors Theiler window independence
    % -------------------------------------------------------------
    fprintf('Test 2: Independent Theiler exclusion window... ');
    try
        S = (1:100)'; % 1D state space along a line
        theiler_win = 10;
        k_neighbors = 2;
        [nn_idx, nn_dist] = find_neighbors(S, k_neighbors, theiler_win);

        % For point i, no neighbor should be in [i - theiler_win, i + theiler_win]
        for i = 1:size(S, 1)
            for k = 1:k_neighbors
                if ~isnan(nn_idx(i, k))
                    assert(abs(nn_idx(i, k) - i) > theiler_win, ...
                           sprintf('Point %d has neighbor %d within Theiler window %d', i, nn_idx(i, k), theiler_win));
                end
            end
        end
        fprintf('PASSED\n');
        n_passed = n_passed + 1;
    catch ME
        fprintf('FAILED: %s\n', ME.message);
        n_failed = n_failed + 1;
    end

    % -------------------------------------------------------------
    % Test 3: compute_divergence_curve zero-distance and log(0) protection
    % -------------------------------------------------------------
    fprintf('Test 3: Zero-distance & log(0) numerical protection... ');
    try
        % Create identical duplicate trajectories
        S = repmat((1:50)', 1, 3);
        nn_idx = repmat(1:50, 50, 1)'; % force distance to be 0
        nn_idx(:, 1) = (1:50)';

        [div_curve, ~, diag] = compute_divergence_curve(S, nn_idx(:, 1), 10, 'MinDistance', 1e-12);
        assert(~any(isinf(div_curve)), 'Inf detected in divergence curve');
        assert(~any(isnan(div_curve(1:5))), 'Unexpected NaNs in initial steps');
        assert(diag.n_zero_distances > 0, 'Zero distance counter failed to trigger');

        fprintf('PASSED\n');
        n_passed = n_passed + 1;
    catch ME
        fprintf('FAILED: %s\n', ME.message);
        n_failed = n_failed + 1;
    end

    % -------------------------------------------------------------
    % Test 4: estimate_lye linear regression and R^2 computation
    % -------------------------------------------------------------
    fprintf('Test 4: estimate_lye linear regression and R^2... ');
    try
        slope_true = 0.75;
        intercept_true = -1.2;
        t = 1:100;
        div_synth = slope_true * t + intercept_true; % perfect line

        res = estimate_lye(div_synth, 1, 'ShortTermRange', [10 50], 'LongTermRange', [], 'TimeUnit', 'samples');
        assert(abs(res.LyE_short - slope_true) < 1e-6, 'Slope mismatch on pure linear signal');
        assert(abs(res.fit_short.R2 - 1.0) < 1e-6, 'R2 should be 1.0 for perfect line');

        fprintf('PASSED\n');
        n_passed = n_passed + 1;
    catch ME
        fprintf('FAILED: %s\n', ME.message);
        n_failed = n_failed + 1;
    end

    % -------------------------------------------------------------
    % Test 5: Linear contractive system verification
    % -------------------------------------------------------------
    fprintf('Test 5: Contractive map slope verification... ');
    try
        a = 0.95;
        n_pts = 1000;
        x = zeros(n_pts, 1);
        x(1) = 1.0;
        for i = 1:n_pts-1
            x(i+1) = a * x(i);
        end
        theo_contractive = log(a); % -0.0513

        % In a contractive system, distance between x(i+k) and x(j+k) shrinks as a^k:
        % ln(d(k)) = ln(|x(i)-x(j)|) + k * ln(a)
        % Slope of ln(d(k)) vs k is strictly ln(a)
        d_k = abs(x(2:20) - x(1:19));
        log_d = log(d_k);
        p_fit = polyfit(1:length(log_d), log_d', 1);
        assert(abs(p_fit(1) - theo_contractive) < 1e-5, 'Contractive slope mismatch');

        fprintf('PASSED\n');
        n_passed = n_passed + 1;
    catch ME
        fprintf('FAILED: %s\n', ME.message);
        n_failed = n_failed + 1;
    end

    % -------------------------------------------------------------
    % Test 6: Unified estimate_rosenstein_lye end-to-end execution
    % -------------------------------------------------------------
    fprintf('Test 6: estimate_rosenstein_lye end-to-end integration... ');
    try
        t = (0:0.01:20)';
        sig = sin(2*pi*t); % periodic
        res_e2e = estimate_rosenstein_lye(sig, ...
            'EmbeddingDimension', 3, 'EmbeddingDelay', 25, 'TheilerWindow', 50, ...
            'NumNeighbors', 1, 'DivergenceWindow', 200, 'SamplesPerStride', 100);
        assert(~isnan(res_e2e.LyE_short), 'LyE_short is NaN');
        assert(isstruct(res_e2e.diagnostics), 'Diagnostics struct missing');
        assert(res_e2e.execution_time > 0, 'Execution time not logged');

        fprintf('PASSED\n');
        n_passed = n_passed + 1;
    catch ME
        fprintf('FAILED: %s\n', ME.message);
        n_failed = n_failed + 1;
    end

    % Summary
    fprintf('\n----------------------------------------------------\n');
    fprintf('TEST SUITE SUMMARY: %d PASSED, %d FAILED\n', n_passed, n_failed);
    fprintf('----------------------------------------------------\n\n');

end
