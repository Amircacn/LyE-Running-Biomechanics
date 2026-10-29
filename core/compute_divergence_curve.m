function [div_curve, div_matrix, diagnostics] = compute_divergence_curve(S, nn_idx, max_steps, varargin)
%COMPUTE_DIVERGENCE_CURVE Track divergence of nearest neighbor pairs.
%
%   [div_curve, div_matrix, diagnostics] = COMPUTE_DIVERGENCE_CURVE(S, ...
%       nn_idx, max_steps) tracks how nearest-neighbor pairs diverge over
%   time in the reconstructed state space. This is the core computation for
%   estimating Lyapunov exponents via the Rosenstein (1993) algorithm.
%
%   [div_curve, div_matrix, diagnostics] = COMPUTE_DIVERGENCE_CURVE(S, ...
%       nn_idx, max_steps, 'Name', Value, ...) specifies additional
%   options as name-value pairs.
%
%   INPUTS
%   ------
%   S          : N x dim state space matrix (from embed_state_space).
%   nn_idx     : N x K matrix of nearest-neighbor indices (from
%                find_neighbors). K is the number of neighbors per point.
%   max_steps  : Maximum number of time steps (positive integer) to track
%                divergence forward in time. The divergence curve will have
%                max_steps elements corresponding to steps 1..max_steps.
%
%   OPTIONAL NAME-VALUE PAIRS
%   -------------------------
%   'MinDistance' : scalar (default: eps)
%       Minimum distance value. Zero distances are replaced with this
%       value before taking the logarithm to avoid log(0) = -Inf.
%
%   OUTPUTS
%   -------
%   div_curve   : 1 x max_steps vector of mean log divergence.
%                 div_curve(k) = mean over valid pairs of
%                 log( ||S(i+k,:) - S(nn_idx(i,j)+k,:)|| )
%                 for time step k = 1, 2, ..., max_steps.
%   div_matrix  : (N*K) x max_steps matrix of raw Euclidean distances
%                 (before log or averaging). Row (i-1)*K + j corresponds
%                 to pair (i, nn_idx(i,j)). Entries are NaN where the
%                 pair cannot be tracked (index out of bounds or NaN
%                 neighbor).
%   diagnostics : struct with quality diagnostics:
%       .n_valid_pairs    : 1 x max_steps vector, number of valid pairs
%                           contributing to the mean at each step.
%       .n_zero_distances : scalar, total count of zero distances
%                           replaced with MinDistance.
%       .n_nan_entries    : scalar, total NaN entries in div_matrix.
%       .n_inf_entries    : scalar, total Inf entries in log(div_matrix).
%
%   ALGORITHM
%   ---------
%   For each reference point i and its neighbor nn_idx(i,j):
%     1. At time step k, compute:
%        d_ij(k) = ||S(i+k, :) - S(nn_idx(i,j)+k, :)||   (Euclidean)
%     2. If d_ij(k) == 0, replace with MinDistance.
%     3. Compute log(d_ij(k)).
%     4. Average log(d_ij(k)) over all valid (i,j) pairs to get div_curve(k).
%   A pair (i,j) is valid at step k if both i+k and nn_idx(i,j)+k are
%   within bounds [1, N] and nn_idx(i,j) is not NaN.
%
%   EXAMPLE
%   -------
%     [S, ~] = embed_state_space(x, 5, 10);
%     [nn_idx, ~] = find_neighbors(S, 1, 50);
%     [div_curve, div_matrix, diag] = compute_divergence_curve(S, nn_idx, 500);
%     plot(div_curve); xlabel('Time step'); ylabel('Mean log divergence');
%
%   REFERENCES
%   ----------
%   Rosenstein MT, Collins JJ, De Luca CJ (1993). A practical method for
%       calculating largest Lyapunov exponents from small data sets.
%       Physica D, 65(1-2), 117-134.
%
%   NOTES
%   -----
%   - Octave-compatible: no implicit expansion used.
%   - The divergence curve is the fundamental input for estimate_lye.
%
%   Version : 1.0
%   Date    : 2024-01-15
%   Author  : LyE Framework

    % ---- Input parsing ----
    p = inputParser;
    p.FunctionName = 'compute_divergence_curve';
    addRequired(p, 'S', @(v) validateattributes(v, {'numeric'}, {'2d', 'nonempty', 'finite'}));
    addRequired(p, 'nn_idx', @(v) validateattributes(v, {'numeric'}, {'2d', 'nonempty'}));
    addRequired(p, 'max_steps', @(v) validateattributes(v, {'numeric'}, {'scalar', 'positive', 'integer'}));
    addParameter(p, 'MinDistance', eps, ...
        @(v) validateattributes(v, {'numeric'}, {'scalar', 'positive', 'finite'}));

    parse(p, S, nn_idx, max_steps, varargin{:});
    min_dist = p.Results.MinDistance;

    [N, dim] = size(S);
    [N_nn, K] = size(nn_idx);

    % ---- Validate compatibility ----
    if N_nn ~= N
        error('compute_divergence_curve:sizeMismatch', ...
              ['nn_idx has %d rows but S has %d rows. They must match ' ...
               '(one neighbor set per state vector).'], N_nn, N);
    end

    if max_steps >= N
        warning('compute_divergence_curve:maxStepsTooLarge', ...
                ['max_steps (%d) >= number of state vectors (%d). ' ...
                 'Later steps will have very few valid pairs.'], max_steps, N);
    end

    % ---- Pre-allocate ----
    n_pairs    = N * K;
    div_matrix = NaN(n_pairs, max_steps);

    % Diagnostic counters
    total_zero = 0;

    % ---- Compute raw distance matrix ----
    % For each pair (i, nn_idx(i,j)), track distance at each step k.
    % Flatten pairs: row index = (i-1)*K + j

    for k = 1:max_steps
        for j = 1:K
            for i = 1:N
                pair_row = (i - 1) * K + j;
                neighbor = nn_idx(i, j);

                % Skip if neighbor is NaN (no valid neighbor found)
                if isnan(neighbor)
                    continue;
                end

                idx_ref  = i + k;
                idx_neig = neighbor + k;

                % Check bounds
                if idx_ref > N || idx_neig > N || idx_ref < 1 || idx_neig < 1
                    continue;
                end

                % Euclidean distance
                d = sqrt(sum((S(idx_ref, :) - S(idx_neig, :)) .^ 2));

                div_matrix(pair_row, k) = d;
            end
        end
    end

    % ---- Handle zero distances ----
    zero_mask = (div_matrix == 0);
    total_zero = sum(zero_mask(:));
    div_matrix(zero_mask) = min_dist;

    % ---- Compute log divergence matrix ----
    log_div_matrix = log(div_matrix);

    % ---- Compute mean log divergence curve ----
    div_curve   = NaN(1, max_steps);
    n_valid_vec = zeros(1, max_steps);

    for k = 1:max_steps
        col = log_div_matrix(:, k);
        valid = isfinite(col);
        n_valid_vec(k) = sum(valid);
        if n_valid_vec(k) > 0
            div_curve(k) = mean(col(valid));
        end
    end

    % ---- Pack diagnostics ----
    diagnostics.n_valid_pairs    = n_valid_vec;
    diagnostics.n_zero_distances = total_zero;
    diagnostics.n_nan_entries    = sum(isnan(div_matrix(:)));
    diagnostics.n_inf_entries    = sum(isinf(log_div_matrix(:)));

end
