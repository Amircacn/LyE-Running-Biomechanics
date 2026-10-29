function [results] = estimate_rosenstein_lye(data, varargin)
%ESTIMATE_ROSENSTEIN_LYE Unified, modular framework for estimating Largest
% Lyapunov Exponent (LLE / LyE) from time series or reconstructed state spaces.
%
%   results = ESTIMATE_ROSENSTEIN_LYE(data) estimates short- and long-term
%   Lyapunov exponents using the revised Rosenstein algorithm with configurable
%   embedding, independent Theiler window, multi-neighbor tracking, and
%   robust diagnostics.
%
%   results = ESTIMATE_ROSENSTEIN_LYE(data, 'Parameter', Value, ...)
%
%   INPUTS
%   ------
%   data : Numeric vector (1D time series) OR N x d matrix (pre-embedded
%          state space). If 1D, state space embedding is performed internally
%          using EmbeddingDimension and EmbeddingDelay. If 2D, the provided
%          state space is used directly.
%
%   OPTIONAL NAME-VALUE PAIRS
%   -------------------------
%   'EmbeddingDimension' : Positive integer (default: 5).
%   'EmbeddingDelay'     : Positive integer in samples (default: 10).
%   'TheilerWindow'      : Non-negative integer in samples (default: 50).
%                          Temporal exclusion window: pairs with |i - j| <= TheilerWindow
%                          are excluded from nearest neighbor selection.
%                          CRUCIAL: TheilerWindow is strictly independent of EmbeddingDelay.
%   'NumNeighbors'       : Positive integer (default: 1).
%                          1 = Classic Rosenstein (1993) algorithm.
%                          >1 = Mehdizadeh (2019) multi-neighbor approach (e.g. 15).
%   'DivergenceWindow'   : Positive integer in samples (default: 1000).
%                          Number of forward time steps to track neighbor divergence.
%   'SamplesPerStride'   : Positive scalar (default: 100).
%                          Number of samples per movement cycle / stride in data.
%   'ShortTermRange'     : 1x2 numeric vector in strides (default: [0 0.5]).
%   'LongTermRange'      : 1x2 numeric vector in strides (default: [4 10]).
%   'TimeUnit'           : 'strides' (default), 'seconds', 'samples', 'model_time'.
%   'FsRaw'              : Sampling frequency in Hz (default: NaN).
%   'ModelDt'            : Step size for ODE systems (default: NaN).
%   'DistanceMetric'     : 'euclidean' (default) or 'max' (Chebyshev).
%   'Normalize'          : Logical (default: false). Whether to z-score normalize 1D data.
%   'MinDistance'        : Positive scalar (default: eps). Floor for zero distances.
%   'PlotResults'        : Logical (default: false).
%   'SavePlot'           : String (default: ''). File path to save diagnostic figure.
%
%   OUTPUTS
%   -------
%   results : Struct containing:
%       .LyE_short       : Estimated short-term Lyapunov exponent.
%       .LyE_long        : Estimated long-term Lyapunov exponent.
%       .LyE_short_unit  : String denoting unit of short-term LyE.
%       .LyE_long_unit   : String denoting unit of long-term LyE.
%       .fit_short       : Fit metrics (slope, intercept, R2, residuals, range_samples).
%       .fit_long        : Fit metrics for long-term region.
%       .div_curve       : 1 x DivergenceWindow mean logarithmic divergence curve.
%       .time_axis       : Time axis values corresponding to div_curve.
%       .diagnostics     : Struct of numerical diagnostics (valid pairs, zeros, NaNs).
%       .params          : Struct recording all operational parameters used.
%       .execution_time  : Total wall-clock time in seconds.
%       .warnings        : Cell array of warnings or potential anomalies encountered.
%
%   REFERENCES
%   ----------
%   - Rosenstein, M. T., Collins, J. J., & De Luca, C. J. (1993). A practical
%     method for calculating largest Lyapunov exponents from small data sets.
%     Physica D, 65(1-2), 117-134.
%   - Mehdizadeh, S. (2019). A robust method to estimate the largest Lyapunov
%     exponent of noisy signals: A revision to the Rosenstein's algorithm.
%     Journal of Biomechanics, 85, 84-91.
%   - Mehdizadeh, S., Arshi, A. R., & Davids, K. (2014). Effect of speed on
%     local dynamic stability of locomotion under different task constraints
%     in running. European Journal of Sport Science, 14(8), 791-798.
%
%   Author: Computational Biomechanics & Dynamical Systems Framework
%   Version: 2.0 (2024)

    t_start = tic;

    % Parse inputs
    p = inputParser;
    p.FunctionName = 'estimate_rosenstein_lye';
    addRequired(p, 'data', @(v) isnumeric(v) && ~isempty(v));
    addParameter(p, 'EmbeddingDimension', 5, @(v) isnumeric(v) && isscalar(v) && v >= 1);
    addParameter(p, 'EmbeddingDelay', 10, @(v) isnumeric(v) && isscalar(v) && v >= 1);
    addParameter(p, 'TheilerWindow', 50, @(v) isnumeric(v) && isscalar(v) && v >= 0);
    addParameter(p, 'NumNeighbors', 1, @(v) isnumeric(v) && isscalar(v) && v >= 1);
    addParameter(p, 'DivergenceWindow', 1000, @(v) isnumeric(v) && isscalar(v) && v >= 1);
    addParameter(p, 'SamplesPerStride', 100, @(v) isnumeric(v) && isscalar(v) && v > 0);
    addParameter(p, 'ShortTermRange', [0 0.5], @(v) isnumeric(v) && (isempty(v) || numel(v) == 2));
    addParameter(p, 'LongTermRange', [4 10], @(v) isnumeric(v) && (isempty(v) || numel(v) == 2));
    addParameter(p, 'TimeUnit', 'strides', @(v) any(strcmpi(v, {'strides', 'seconds', 'samples', 'model_time'})));
    addParameter(p, 'FsRaw', NaN, @(v) isnumeric(v) && isscalar(v));
    addParameter(p, 'ModelDt', NaN, @(v) isnumeric(v) && isscalar(v));
    addParameter(p, 'DistanceMetric', 'euclidean', @(v) any(strcmpi(v, {'euclidean', 'max'})));
    addParameter(p, 'Normalize', false, @(v) islogical(v) || isnumeric(v));
    addParameter(p, 'MinDistance', eps, @(v) isnumeric(v) && isscalar(v) && v > 0);
    addParameter(p, 'PlotResults', false, @(v) islogical(v) || isnumeric(v));
    addParameter(p, 'SavePlot', '', @ischar);

    parse(p, data, varargin{:});
    opts = p.Results;

    warn_list = {};

    % Step 1: Reconstruct or Validate State Space
    if isvector(data)
        x = data(:);
        [S, embed_params] = embed_state_space(x, opts.EmbeddingDimension, opts.EmbeddingDelay, ...
                                              'Normalize', opts.Normalize);
    else
        S = data;
        embed_params.dim = size(S, 2);
        embed_params.delay = NaN;
        embed_params.normalize = false;
        embed_params.input_length = size(S, 1);
        embed_params.n_vectors = size(S, 1);
        embed_params.mean = NaN;
        embed_params.std = NaN;
    end

    [N_vectors, dim] = size(S);

    % Sanity check on divergence window
    actual_max_steps = opts.DivergenceWindow;
    if actual_max_steps >= N_vectors
        warn_msg = sprintf('DivergenceWindow (%d) >= available state vectors (%d). Truncating to %d.', ...
                           actual_max_steps, N_vectors, N_vectors - 1);
        warn_list{end+1} = warn_msg;
        actual_max_steps = max(1, N_vectors - 1);
    end

    % Step 2: Neighbor Identification with Independent Theiler Window
    [nn_idx, nn_dist] = find_neighbors(S, opts.NumNeighbors, opts.TheilerWindow, ...
                                       'DistanceMetric', opts.DistanceMetric);

    % Check if any points lacked valid neighbors
    nan_neighbor_count = sum(isnan(nn_idx(:)));
    if nan_neighbor_count > 0
        warn_msg = sprintf('%d neighbor slots could not be populated (Theiler window too large or boundary).', nan_neighbor_count);
        warn_list{end+1} = warn_msg;
    end

    % Step 3: Divergence Tracking
    [div_curve, ~, diagnostics] = compute_divergence_curve(S, nn_idx, actual_max_steps, ...
                                                           'MinDistance', opts.MinDistance);

    % Step 4: LyE Estimation (linear regression over defined windows)
    est_results = estimate_lye(div_curve, opts.SamplesPerStride, ...
                               'ShortTermRange', opts.ShortTermRange, ...
                               'LongTermRange', opts.LongTermRange, ...
                               'FsRaw', opts.FsRaw, ...
                               'TimeUnit', opts.TimeUnit, ...
                               'PlotResults', opts.PlotResults, ...
                               'SavePlot', opts.SavePlot);

    % Consolidate warnings
    if ~isempty(est_results.warnings)
        warn_list = [warn_list, est_results.warnings];
    end

    % Pack final results struct
    results = est_results;
    results.diagnostics = diagnostics;
    results.params = opts;
    results.params.actual_max_steps = actual_max_steps;
    results.params.embedding_info = embed_params;
    results.execution_time = toc(t_start);
    results.warnings = warn_list;

end
