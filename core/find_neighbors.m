function [nn_idx, nn_dist] = find_neighbors(S, num_neighbors, theiler_window, varargin)
%FIND_NEIGHBORS Find nearest neighbors in state space with Theiler window.
%
%   [nn_idx, nn_dist] = FIND_NEIGHBORS(S, num_neighbors, theiler_window)
%   finds the num_neighbors nearest neighbors for each point in the state
%   space matrix S, excluding temporally close points within a Theiler
%   window. The Theiler window is INDEPENDENT from the embedding delay.
%
%   [nn_idx, nn_dist] = FIND_NEIGHBORS(S, num_neighbors, theiler_window, ...
%                                       'Name', Value, ...)
%   specifies additional options as name-value pairs.
%
%   INPUTS
%   ------
%   S               : N x dim state space matrix (each row is a state
%                     vector, from embed_state_space).
%   num_neighbors   : Number of nearest neighbors to find per point
%                     (positive integer). Use 1 for classic Rosenstein
%                     algorithm, >1 for Mehdizadeh-type methods.
%   theiler_window  : Temporal separation window in samples (non-negative
%                     integer). Points within [i - theiler_window,
%                     i + theiler_window] are excluded from being neighbors
%                     of point i. This prevents selecting temporally
%                     correlated points as neighbors.
%                     NOTE: This parameter is INDEPENDENT of the embedding
%                     delay (tau). They serve different purposes.
%
%   OPTIONAL NAME-VALUE PAIRS
%   -------------------------
%   'DistanceMetric' : string (default: 'euclidean')
%       Distance metric for neighbor search:
%         'euclidean' - L2 norm (standard Euclidean distance)
%         'max'       - Chebyshev / L-infinity norm (max absolute
%                       difference across dimensions)
%
%   OUTPUTS
%   -------
%   nn_idx  : N x num_neighbors matrix of neighbor indices. nn_idx(i,k) is
%             the index in S of the k-th nearest neighbor of point i.
%             If fewer than num_neighbors valid neighbors exist for point i,
%             the remaining entries are NaN.
%   nn_dist : N x num_neighbors matrix of corresponding distances.
%             Entries corresponding to NaN indices are also NaN.
%
%   EXAMPLE
%   -------
%     % Find 1 nearest neighbor with Theiler window of 50 samples
%     [S, ~] = embed_state_space(x, 5, 10);
%     [nn_idx, nn_dist] = find_neighbors(S, 1, 50);
%
%     % Find 5 nearest neighbors using Chebyshev distance
%     [nn_idx, nn_dist] = find_neighbors(S, 5, 50, 'DistanceMetric', 'max');
%
%   REFERENCES
%   ----------
%   Rosenstein MT, Collins JJ, De Luca CJ (1993). A practical method for
%       calculating largest Lyapunov exponents from small data sets.
%       Physica D, 65(1-2), 117-134.
%   Theiler J (1986). Spurious dimension from correlation algorithms
%       applied to limited time-series data. Physical Review A, 34(3).
%
%   NOTES
%   -----
%   - For large state spaces this function processes point-by-point to
%     manage memory. For moderate sizes (N < ~10000) it uses a full
%     distance matrix for speed.
%   - Octave-compatible: uses bsxfun instead of implicit expansion.
%
%   Version : 1.0
%   Date    : 2024-01-15
%   Author  : LyE Framework

    % ---- Input parsing ----
    p = inputParser;
    p.FunctionName = 'find_neighbors';
    addRequired(p, 'S', @(v) validateattributes(v, {'numeric'}, {'2d', 'nonempty', 'finite'}));
    addRequired(p, 'num_neighbors', @(v) validateattributes(v, {'numeric'}, {'scalar', 'positive', 'integer'}));
    addRequired(p, 'theiler_window', @(v) validateattributes(v, {'numeric'}, {'scalar', 'nonnegative', 'integer'}));
    addParameter(p, 'DistanceMetric', 'euclidean', ...
        @(v) any(strcmpi(v, {'euclidean', 'max'})));

    parse(p, S, num_neighbors, theiler_window, varargin{:});
    metric = lower(p.Results.DistanceMetric);

    [N, dim] = size(S);

    if N < 2
        error('find_neighbors:tooFewPoints', ...
              'State space must have at least 2 points, got %d.', N);
    end

    % ---- Pre-allocate outputs ----
    nn_idx  = NaN(N, num_neighbors);
    nn_dist = NaN(N, num_neighbors);

    % ---- Strategy selection based on N ----
    % For moderate N, compute full pairwise distance matrix (fast).
    % For large N, iterate point-by-point (memory efficient).
    use_full_matrix = (N <= 10000);

    if use_full_matrix
        % Compute full N x N distance matrix
        D = compute_distance_matrix(S, metric);

        for i = 1:N
            % Zero out (set to Inf) points within Theiler window
            dists = D(i, :);
            excl_lo = max(1, i - theiler_window);
            excl_hi = min(N, i + theiler_window);
            dists(excl_lo:excl_hi) = Inf;

            % Find k smallest distances
            [sorted_d, sorted_i] = sort(dists, 'ascend');

            % How many valid neighbors are available?
            n_valid = sum(isfinite(sorted_d));
            n_use   = min(num_neighbors, n_valid);

            if n_use > 0
                nn_idx(i, 1:n_use)  = sorted_i(1:n_use);
                nn_dist(i, 1:n_use) = sorted_d(1:n_use);
            end
        end
    else
        % Point-by-point computation for large N
        for i = 1:N
            % Compute distances from point i to all other points
            % diff: N x dim
            diff_mat = bsxfun(@minus, S, S(i, :));

            if strcmp(metric, 'euclidean')
                dists = sqrt(sum(diff_mat .^ 2, 2))';
            else  % 'max' (Chebyshev)
                dists = max(abs(diff_mat), [], 2)';
            end

            % Exclude Theiler window
            excl_lo = max(1, i - theiler_window);
            excl_hi = min(N, i + theiler_window);
            dists(excl_lo:excl_hi) = Inf;

            % Find k smallest
            [sorted_d, sorted_i] = sort(dists, 'ascend');

            n_valid = sum(isfinite(sorted_d));
            n_use   = min(num_neighbors, n_valid);

            if n_use > 0
                nn_idx(i, 1:n_use)  = sorted_i(1:n_use);
                nn_dist(i, 1:n_use) = sorted_d(1:n_use);
            end
        end
    end

end

% =========================================================================
% LOCAL FUNCTION: compute full pairwise distance matrix
% =========================================================================
function D = compute_distance_matrix(S, metric)
%COMPUTE_DISTANCE_MATRIX Compute N x N pairwise distance matrix.
    [N, ~] = size(S);

    if strcmp(metric, 'euclidean')
        % Efficient Euclidean: ||a-b||^2 = ||a||^2 + ||b||^2 - 2*a.b
        sq_norms = sum(S .^ 2, 2);  % N x 1
        % D_sq = sq_norms + sq_norms' - 2*(S*S')
        D_sq = bsxfun(@plus, sq_norms, sq_norms') - 2 * (S * S');
        % Clamp small negative values from floating-point errors
        D_sq(D_sq < 0) = 0;
        D = sqrt(D_sq);
    else  % 'max' (Chebyshev)
        D = zeros(N, N);
        for i = 1:N
            diff_mat = bsxfun(@minus, S, S(i, :));
            D(i, :) = max(abs(diff_mat), [], 2)';
        end
    end
end
