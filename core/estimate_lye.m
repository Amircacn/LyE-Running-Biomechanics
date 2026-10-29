function [results] = estimate_lye(div_curve, samples_per_stride, varargin)
%ESTIMATE_LYE Estimate Lyapunov exponents from a divergence curve.
%
%   results = ESTIMATE_LYE(div_curve, samples_per_stride) estimates
%   short-term and long-term Lyapunov exponents by fitting linear
%   regressions to specified regions of the mean log divergence curve.
%
%   results = ESTIMATE_LYE(div_curve, samples_per_stride, 'Name', Value, ...)
%   specifies additional options as name-value pairs.
%
%   INPUTS
%   ------
%   div_curve          : 1 x T vector of mean log divergence values
%                        (from compute_divergence_curve). T is the number
%                        of time steps tracked.
%   samples_per_stride : Number of samples per stride/cycle in normalized
%                        data (positive scalar). For gait data normalized
%                        to 100 samples per stride, use 100. For non-gait
%                        data, this defines the time unit conversion factor
%                        (samples per one "unit").
%
%   OPTIONAL NAME-VALUE PAIRS
%   -------------------------
%   'ShortTermRange' : [start_stride end_stride] (default: [0 0.5])
%       Range in strides for short-term LyE fitting. Converted to sample
%       indices via floor/ceil(value * samples_per_stride).
%       Start is inclusive, end is inclusive (in strides).
%
%   'LongTermRange'  : [start_stride end_stride] (default: [4 10])
%       Range in strides for long-term LyE fitting. Same conversion.
%
%   'FsRaw'          : scalar (default: [])
%       Raw sampling frequency in Hz. If provided along with
%       TimeUnit='seconds', the time axis is converted to physical time.
%       Not used unless TimeUnit is 'seconds'.
%
%   'TimeUnit'       : string (default: 'strides')
%       Unit for the output time axis and exponent units:
%         'strides'  - time in strides (sample_index / samples_per_stride)
%         'samples'  - time in samples (raw index)
%         'seconds'  - time in seconds (requires FsRaw)
%
%   'PlotResults'    : logical (default: false)
%       If true, generate a diagnostic plot showing the divergence curve
%       and fitted regression lines.
%
%   'SavePlot'       : string (default: '')
%       File path to save the plot (e.g., 'results/lye_plot.png'). If
%       empty, the plot is not saved. Ignored if PlotResults is false.
%
%   OUTPUTS
%   -------
%   results : struct with fields:
%     .LyE_short      : Short-term Lyapunov exponent (slope of linear fit)
%     .LyE_long       : Long-term Lyapunov exponent
%     .LyE_short_unit : String describing units (e.g., 'per stride')
%     .LyE_long_unit  : String describing units
%     .fit_short       : struct with:
%         .slope          - slope of linear regression
%         .intercept      - y-intercept
%         .R2             - coefficient of determination
%         .residuals      - residual vector
%         .range_samples  - [start_sample end_sample] indices used
%     .fit_long        : same structure as fit_short
%     .div_curve       : the input divergence curve
%     .time_axis       : time axis in the chosen unit (1 x T)
%     .warnings        : cell array of warning strings generated
%
%   ALGORITHM
%   ---------
%   1. Convert the stride-based fitting ranges to sample indices.
%   2. Extract the portion of div_curve within each range.
%   3. Fit a least-squares linear regression (polyfit degree 1).
%   4. The slope of the line is the Lyapunov exponent estimate.
%   5. Compute R^2 and residuals for fit quality assessment.
%
%   The Lyapunov exponent is the slope in the chosen time unit:
%     - 'strides': slope of log(divergence) vs time-in-strides
%     - 'samples': slope of log(divergence) vs time-in-samples
%     - 'seconds': slope of log(divergence) vs time-in-seconds
%
%   EXAMPLE
%   -------
%     % Estimate LyE from a divergence curve (gait data, 100 samples/stride)
%     results = estimate_lye(div_curve, 100);
%     fprintf('Short-term LyE: %.4f %s\n', results.LyE_short, results.LyE_short_unit);
%
%     % Custom ranges with plotting
%     results = estimate_lye(div_curve, 100, ...
%         'ShortTermRange', [0 1], 'LongTermRange', [4 8], ...
%         'PlotResults', true);
%
%     % Time axis in seconds (raw sampling at 120 Hz)
%     results = estimate_lye(div_curve, 100, ...
%         'FsRaw', 120, 'TimeUnit', 'seconds', 'PlotResults', true);
%
%   REFERENCES
%   ----------
%   Rosenstein MT, Collins JJ, De Luca CJ (1993). A practical method for
%       calculating largest Lyapunov exponents from small data sets.
%       Physica D, 65(1-2), 117-134.
%   Mehdizadeh S (2018). The largest Lyapunov exponent of gait in young
%       and older adults: A systematic review and meta-analysis.
%       Gait & Posture, 60, 241-250.
%
%   NOTES
%   -----
%   - Octave-compatible.
%   - If a fitting range exceeds available data, the exponent is NaN and
%     a warning is recorded.
%
%   Version : 1.0
%   Date    : 2024-01-15
%   Author  : LyE Framework

    % ---- Input parsing ----
    p = inputParser;
    p.FunctionName = 'estimate_lye';
    addRequired(p, 'div_curve', @(v) validateattributes(v, {'numeric'}, {'vector', 'nonempty'}));
    addRequired(p, 'samples_per_stride', @(v) validateattributes(v, {'numeric'}, {'scalar', 'positive', 'finite'}));
    addParameter(p, 'ShortTermRange', [0 0.5], @(v) validate_range(v));
    addParameter(p, 'LongTermRange',  [4 10],  @(v) validate_range(v));
    addParameter(p, 'FsRaw', [], @(v) isempty(v) || (isnumeric(v) && isscalar(v) && (isnan(v) || v > 0)));
    addParameter(p, 'TimeUnit', 'strides', ...
        @(v) any(strcmpi(v, {'strides', 'samples', 'seconds'})));
    addParameter(p, 'PlotResults', false, @(v) islogical(v) || isnumeric(v));
    addParameter(p, 'SavePlot', '', @ischar);

    parse(p, div_curve, samples_per_stride, varargin{:});

    short_range  = p.Results.ShortTermRange;
    long_range   = p.Results.LongTermRange;
    fs_raw       = p.Results.FsRaw;
    time_unit    = lower(p.Results.TimeUnit);
    do_plot      = logical(p.Results.PlotResults);
    save_path    = p.Results.SavePlot;

    % Force row vector
    div_curve = div_curve(:)';
    T = length(div_curve);

    % Collect warnings
    warn_list = {};

    % ---- Validate TimeUnit = 'seconds' requires FsRaw ----
    if strcmp(time_unit, 'seconds') && isempty(fs_raw)
        error('estimate_lye:noFsRaw', ...
              'TimeUnit is ''seconds'' but FsRaw was not provided.');
    end

    % ---- Build time axis ----
    sample_indices = 1:T;  % 1-based sample index

    switch time_unit
        case 'strides'
            time_axis = sample_indices / samples_per_stride;
            unit_str  = 'per stride';
        case 'samples'
            time_axis = sample_indices;
            unit_str  = 'per sample';
        case 'seconds'
            % Convert sample steps to seconds
            % Each sample step corresponds to 1/fs_raw seconds in the
            % original data (before any normalization)
            % But since we have samples_per_stride for stride normalization,
            % if the data was resampled to samples_per_stride, the
            % effective dt depends on context.
            % For normalized data: stride_duration = samples_per_stride / fs_raw is NOT right
            % Instead: time_in_seconds = sample_index / fs_raw
            % But sample_index here is after embedding, so each step is
            % one sample in the resampled series.
            % If FsRaw is the frequency of the embedded time series:
            time_axis = sample_indices / fs_raw;
            unit_str  = 'per second';
    end

    % ---- Convert stride ranges to sample indices & Fit ----
    if ~isempty(short_range)
        short_start_samp = max(1, floor(short_range(1) * samples_per_stride) + 1);
        short_end_samp   = ceil(short_range(2) * samples_per_stride);
        short_start_samp = max(1, short_start_samp);
        [fit_short, w] = fit_region(div_curve, time_axis, short_start_samp, short_end_samp, T, 'Short-term');
        warn_list = [warn_list, w];
    else
        fit_short.slope = NaN;
        fit_short.intercept = NaN;
        fit_short.R2 = NaN;
        fit_short.residuals = [];
        fit_short.range_samples = [];
    end

    if ~isempty(long_range)
        long_start_samp  = max(1, floor(long_range(1) * samples_per_stride) + 1);
        long_end_samp    = ceil(long_range(2) * samples_per_stride);
        long_start_samp  = max(1, long_start_samp);
        [fit_long, w] = fit_region(div_curve, time_axis, long_start_samp, long_end_samp, T, 'Long-term');
        warn_list = [warn_list, w];
    else
        fit_long.slope = NaN;
        fit_long.intercept = NaN;
        fit_long.R2 = NaN;
        fit_long.residuals = [];
        fit_long.range_samples = [];
    end

    % ---- Pack results ----
    results.LyE_short      = fit_short.slope;
    results.LyE_long       = fit_long.slope;
    results.LyE_short_unit = unit_str;
    results.LyE_long_unit  = unit_str;
    results.fit_short       = fit_short;
    results.fit_long        = fit_long;
    results.div_curve       = div_curve;
    results.time_axis       = time_axis;
    results.warnings        = warn_list;

    % ---- Optional plotting ----
    if do_plot
        plot_divergence(results, time_unit, save_path);
    end

end

% =========================================================================
% LOCAL FUNCTION: fit linear regression to a region
% =========================================================================
function [fit_out, warnings] = fit_region(div_curve, time_axis, start_samp, end_samp, T, label)
%FIT_REGION Fit a linear regression to a region of the divergence curve.

    warnings = {};
    fit_out.slope         = NaN;
    fit_out.intercept     = NaN;
    fit_out.R2            = NaN;
    fit_out.residuals     = [];
    fit_out.range_samples = [start_samp, end_samp];

    % ---- Check if range is valid ----
    if start_samp > T
        warnings{end+1} = sprintf('%s range starts at sample %d but only %d steps available. Returning NaN.', ...
                                  label, start_samp, T);
        return;
    end

    if end_samp > T
        warnings{end+1} = sprintf('%s range ends at sample %d but only %d steps available. Truncating to %d.', ...
                                  label, end_samp, T, T);
        end_samp = T;
        fit_out.range_samples(2) = end_samp;
    end

    if start_samp > end_samp
        warnings{end+1} = sprintf('%s range is empty (start=%d > end=%d). Returning NaN.', ...
                                  label, start_samp, end_samp);
        return;
    end

    % ---- Extract region ----
    idx = start_samp:end_samp;
    y = div_curve(idx);
    t = time_axis(idx);

    % Remove NaN values
    valid = isfinite(y) & isfinite(t);
    y = y(valid);
    t = t(valid);

    if length(y) < 2
        warnings{end+1} = sprintf('%s region has fewer than 2 valid data points. Returning NaN.', label);
        return;
    end

    % ---- Linear regression via polyfit ----
    coeffs = polyfit(t, y, 1);
    fit_out.slope     = coeffs(1);
    fit_out.intercept = coeffs(2);

    % Predicted values
    y_pred = polyval(coeffs, t);
    fit_out.residuals = y - y_pred;

    % R-squared
    SS_res = sum(fit_out.residuals .^ 2);
    SS_tot = sum((y - mean(y)) .^ 2);
    if SS_tot == 0
        fit_out.R2 = NaN;
        warnings{end+1} = sprintf('%s divergence curve is constant in fitting region. R^2 undefined.', label);
    else
        fit_out.R2 = 1 - SS_res / SS_tot;
    end

end

% =========================================================================
% LOCAL FUNCTION: validate range input
% =========================================================================
function tf = validate_range(v)
%VALIDATE_RANGE Check that v is a 2-element numeric vector with v(1)<=v(2) or empty.
    if isempty(v)
        tf = true;
        return;
    end
    tf = isnumeric(v) && numel(v) == 2 && v(1) >= 0 && v(2) >= v(1);
    if ~tf
        error('estimate_lye:invalidRange', ...
              'Range must be a 2-element vector [start end] with 0 <= start <= end, or empty.');
    end
end

% =========================================================================
% LOCAL FUNCTION: plot divergence curve with fits
% =========================================================================
function plot_divergence(results, time_unit, save_path)
%PLOT_DIVERGENCE Create diagnostic plot of divergence curve and fits.

    fig = figure('Name', 'Lyapunov Exponent Estimation', ...
                 'NumberTitle', 'off', 'Color', 'w');

    % Divergence curve
    plot(results.time_axis, results.div_curve, 'b-', 'LineWidth', 1);
    hold on;

    % Short-term fit region
    if ~isnan(results.fit_short.slope)
        rs = results.fit_short.range_samples;
        idx_s = rs(1):rs(2);
        idx_s = idx_s(idx_s <= length(results.time_axis));
        if ~isempty(idx_s)
            t_fit = results.time_axis(idx_s);
            y_fit = results.fit_short.intercept + results.fit_short.slope * t_fit;
            plot(t_fit, y_fit, 'r-', 'LineWidth', 2.5);
        end
    end

    % Long-term fit region
    if ~isnan(results.fit_long.slope)
        rl = results.fit_long.range_samples;
        idx_l = rl(1):rl(2);
        idx_l = idx_l(idx_l <= length(results.time_axis));
        if ~isempty(idx_l)
            t_fit = results.time_axis(idx_l);
            y_fit = results.fit_long.intercept + results.fit_long.slope * t_fit;
            plot(t_fit, y_fit, 'g-', 'LineWidth', 2.5);
        end
    end

    hold off;

    % Labels
    switch time_unit
        case 'strides'
            xlabel('Time (strides)');
        case 'samples'
            xlabel('Time (samples)');
        case 'seconds'
            xlabel('Time (seconds)');
    end
    ylabel('Mean log divergence');
    title('Lyapunov Exponent Estimation');

    % Legend
    leg_entries = {'Divergence curve'};
    if ~isnan(results.fit_short.slope)
        leg_entries{end+1} = sprintf('Short-term LyE = %.4f (R^2=%.3f)', ...
                                     results.LyE_short, results.fit_short.R2);
    end
    if ~isnan(results.fit_long.slope)
        leg_entries{end+1} = sprintf('Long-term LyE = %.4f (R^2=%.3f)', ...
                                     results.LyE_long, results.fit_long.R2);
    end
    legend(leg_entries, 'Location', 'best');

    grid on;
    set(gca, 'FontSize', 11);

    % Save if requested
    if ~isempty(save_path)
        % Ensure directory exists
        [save_dir, ~, ~] = fileparts(save_path);
        if ~isempty(save_dir) && ~exist(save_dir, 'dir')
            mkdir(save_dir);
        end
        saveas(fig, save_path);
        fprintf('Plot saved to: %s\n', save_path);
    end

end
