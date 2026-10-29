function [stride_indices, cycle_info] = detect_gait_cycles(signal, fs, varargin)
%DETECT_GAIT_CYCLES Robust automated gait cycle / stride segmentation from
% continuous biomechanical kinematic or kinetic time series.
%
%   [stride_indices, cycle_info] = DETECT_GAIT_CYCLES(signal, fs)
%   detects stride boundary events (e.g. heelstrikes, knee extension/flexion peaks)
%   using prominence-based peak/valley detection tailored for human locomotion.
%
%   INPUTS
%   ------
%   signal : 1D numeric vector (joint angle, velocity, or marker trajectory).
%   fs     : Sampling frequency in Hz (e.g. 100 or 200 Hz).
%
%   OPTIONAL PARAMETERS
%   -------------------
%   'EventType'      : 'peaks' (default) or 'valleys'.
%   'MinStrideTime'  : Minimum expected stride duration in seconds (default: 0.45 s for running).
%   'MaxStrideTime'  : Maximum expected stride duration in seconds (default: 1.5 s).
%   'ProminenceFrac' : Minimum peak prominence as fraction of signal range (default: 0.15).
%
%   OUTPUTS
%   -------
%   stride_indices : Vector of sample indices where stride events occur.
%   cycle_info     : Struct with duration statistics (mean, std, stride_times).
%
%   Author: Computational Biomechanics & Dynamical Systems Framework
%   Version: 2.0 (2024)

    p = inputParser;
    addRequired(p, 'signal', @(v) isnumeric(v) && isvector(v) && ~isempty(v));
    addRequired(p, 'fs', @(v) isnumeric(v) && isscalar(v) && v > 0);
    addParameter(p, 'EventType', 'peaks', @(v) any(strcmpi(v, {'peaks', 'valleys'})));
    addParameter(p, 'MinStrideTime', 0.45, @(v) isnumeric(v) && isscalar(v) && v > 0);
    addParameter(p, 'MaxStrideTime', 1.5, @(v) isnumeric(v) && isscalar(v) && v > 0);
    addParameter(p, 'ProminenceFrac', 0.15, @(v) isnumeric(v) && isscalar(v) && v > 0);

    parse(p, signal, fs, varargin{:});
    sig = signal(:);
    event_type = lower(p.Results.EventType);
    min_dist_samples = round(p.Results.MinStrideTime * fs);
    min_prom = p.Results.ProminenceFrac * (max(sig) - min(sig));

    if strcmp(event_type, 'valleys')
        [~, locs] = findpeaks(-sig, 'MinPeakDistance', min_dist_samples, 'MinPeakProminence', min_prom);
    else
        [~, locs] = findpeaks(sig, 'MinPeakDistance', min_dist_samples, 'MinPeakProminence', min_prom);
    end

    if length(locs) < 2
        warning('detect_gait_cycles:tooFewStrides', 'Fewer than 2 strides detected.');
        stride_indices = locs;
        cycle_info.n_strides = 0;
        cycle_info.mean_stride_time = NaN;
        cycle_info.std_stride_time = NaN;
        cycle_info.stride_times = [];
        return;
    end

    stride_times = diff(locs) / fs;

    % Filter outlier stride durations
    valid_mask = (stride_times >= p.Results.MinStrideTime) & (stride_times <= p.Results.MaxStrideTime);
    valid_locs = locs([true; valid_mask]);

    stride_indices = valid_locs;
    cycle_info.n_strides = length(stride_indices) - 1;
    cycle_info.stride_times = diff(stride_indices) / fs;
    cycle_info.mean_stride_time = mean(cycle_info.stride_times);
    cycle_info.std_stride_time = std(cycle_info.stride_times);
    cycle_info.cadence_spm = 60 / cycle_info.mean_stride_time * 2; % steps per minute

end
