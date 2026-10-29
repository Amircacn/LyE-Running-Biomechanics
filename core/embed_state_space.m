function [S, params] = embed_state_space(x, dim, delay, varargin)
%EMBED_STATE_SPACE Reconstruct state space using time-delay embedding.
%
%   [S, params] = EMBED_STATE_SPACE(x, dim, delay) reconstructs the state
%   space from a 1-D time series x using time-delay (Takens) embedding with
%   embedding dimension dim and embedding delay (lag) delay.
%
%   [S, params] = EMBED_STATE_SPACE(x, dim, delay, 'Name', Value, ...)
%   specifies additional options as name-value pairs.
%
%   INPUTS
%   ------
%   x      : 1-D numeric time series (row or column vector). Converted to
%             a column vector internally.
%   dim    : Embedding dimension (positive integer >= 1).
%   delay  : Embedding delay in samples (positive integer >= 1). Also known
%             as the embedding lag or tau.
%
%   OPTIONAL NAME-VALUE PAIRS
%   -------------------------
%   'Normalize' : logical (default: false)
%       If true, z-score normalize the time series before embedding.
%       Normalization uses (x - mean(x)) / std(x).
%
%   OUTPUTS
%   -------
%   S      : (N - (dim-1)*delay) x dim matrix of embedded state vectors.
%             Row i of S is [x(i), x(i+delay), x(i+2*delay), ...,
%             x(i+(dim-1)*delay)].
%   params : struct containing all parameters used for reconstruction:
%             .dim          - embedding dimension
%             .delay        - embedding delay (samples)
%             .normalize    - whether normalization was applied
%             .input_length - length of input time series
%             .n_vectors    - number of state vectors produced
%             .mean         - mean of input (NaN if not normalized)
%             .std          - std of input (NaN if not normalized)
%
%   EXAMPLE
%   -------
%     % Embed a sine wave with dim=3, delay=10
%     t = (0:0.01:10)';
%     x = sin(2*pi*t);
%     [S, params] = embed_state_space(x, 3, 10);
%     plot3(S(:,1), S(:,2), S(:,3));
%
%     % With normalization
%     [S, params] = embed_state_space(x, 3, 10, 'Normalize', true);
%
%   REFERENCES
%   ----------
%   Takens F (1981). Detecting strange attractors in turbulence.
%       Lecture Notes in Mathematics, 898, 366-381.
%
%   NOTES
%   -----
%   - Octave-compatible: no MATLAB-only syntax used.
%   - Both row and column vectors are accepted; output S is always
%     (N_vectors x dim).
%
%   Version : 1.0
%   Date    : 2024-01-15
%   Author  : LyE Framework

    % ---- Input parsing ----
    p = inputParser;
    p.FunctionName = 'embed_state_space';
    addRequired(p, 'x', @(v) validateattributes(v, {'numeric'}, {'vector', 'nonempty', 'finite'}));
    addRequired(p, 'dim', @(v) validateattributes(v, {'numeric'}, {'scalar', 'positive', 'integer'}));
    addRequired(p, 'delay', @(v) validateattributes(v, {'numeric'}, {'scalar', 'positive', 'integer'}));
    addParameter(p, 'Normalize', false, @(v) islogical(v) || isnumeric(v));

    parse(p, x, dim, delay, varargin{:});
    do_normalize = logical(p.Results.Normalize);

    % ---- Force column vector ----
    x = x(:);
    N = length(x);

    % ---- Validate sufficient length ----
    min_length = (dim - 1) * delay + 1;
    if N < min_length
        error('embed_state_space:insufficientLength', ...
              ['Input time series has %d samples but embedding requires at least %d ' ...
               'samples (dim=%d, delay=%d).'], N, min_length, dim, delay);
    end

    % ---- Optional normalization ----
    x_mean = NaN;
    x_std  = NaN;
    if do_normalize
        x_mean = mean(x);
        x_std  = std(x);
        if x_std == 0
            error('embed_state_space:zeroVariance', ...
                  'Cannot z-score normalize: time series has zero variance.');
        end
        x = (x - x_mean) / x_std;
    end

    % ---- Time-delay embedding ----
    n_vectors = N - (dim - 1) * delay;

    % Build index matrix: each row selects the indices for one state vector
    % Column j uses indices offset by (j-1)*delay
    col_offsets = (0:dim-1) * delay;          % 1 x dim
    row_starts  = (1:n_vectors)';             % n_vectors x 1
    idx_matrix  = bsxfun(@plus, row_starts, col_offsets);  % n_vectors x dim

    S = x(idx_matrix);

    % ---- Pack output parameters ----
    params.dim          = dim;
    params.delay        = delay;
    params.normalize    = do_normalize;
    params.input_length = N;
    params.n_vectors    = n_vectors;
    params.mean         = x_mean;
    params.std          = x_std;

end
