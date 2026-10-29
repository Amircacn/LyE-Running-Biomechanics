function benchmarks = generate_benchmarks(varargin)
% GENERATE_BENCHMARKS  Generate synthetic benchmark datasets for Lyapunov
% exponent estimation algorithm validation.
%
%   benchmarks = GENERATE_BENCHMARKS()
%       Generates all four benchmark datasets (A-D) with default parameters,
%       saves to benchmarks.mat, and returns the struct.
%
%   benchmarks = GENERATE_BENCHMARKS('Seeds', [42 43 44 45 46])
%       Use custom seeds for independent realizations.
%
%   benchmarks = GENERATE_BENCHMARKS('Lengths', [1000 4000 8000 16000])
%       Use custom data lengths (number of sample points).
%
%   benchmarks = GENERATE_BENCHMARKS('OutputFile', 'my_benchmarks.mat')
%       Save to a custom output file path.
%
% DATASETS:
%   A) Periodic (Van der Pol oscillator, mu=1) — Theoretical LyE = 0
%   B) Lorenz-63 chaotic attractor              — Reference LyE ≈ 0.9056
%   C) Logistic map at r=4                      — Theoretical LyE = ln(2)
%   D) Linear contractive map (a=0.95)          — Theoretical LyE = ln(0.95)
%
% OUTPUT STRUCT FIELDS:
%   benchmarks.periodic     — Van der Pol limit cycle data
%   benchmarks.lorenz       — Lorenz-63 chaotic data (x, y, z components)
%   benchmarks.logistic     — Logistic map data
%   benchmarks.contractive  — Linear contractive map data
%
%   Each sub-struct contains:
%     .clean{realization}{length_idx}   — clean time series
%     .noisy{snr_idx}{realization}{length_idx} — noise-corrupted versions
%     .params        — parameter struct documenting all generation settings
%     .theoretical_lye — theoretical/reference Lyapunov exponent value
%     .time_unit     — string: 'model_time_unit' or 'iteration'
%     .seeds         — array of RNG seeds used
%     .lengths       — array of data lengths generated
%     .snr_levels    — array of SNR levels in dB
%
% REFERENCES:
%   Sprott, J.C. (2003). Chaos and Time-Series Analysis. Oxford Univ. Press.
%   Viswanath, D. (1998). Lyapunov exponents from random Fibonacci
%       sequences to the Lorenz equations. PhD thesis, Cornell University.
%   Strogatz, S.H. (2015). Nonlinear Dynamics and Chaos. Westview Press.
%
% Author: LyPE Benchmark Suite
% Date:   2024

    %% ====================================================================
    %  PARSE INPUT ARGUMENTS
    % =====================================================================
    p = inputParser;
    p.addParameter('Seeds',      [42 43 44 45 46], @(x) isnumeric(x) && isvector(x));
    p.addParameter('Lengths',    [1000 4000 8000 16000], @(x) isnumeric(x) && isvector(x));
    p.addParameter('OutputFile', fullfile(fileparts(mfilename('fullpath')), 'benchmarks.mat'), @ischar);
    p.parse(varargin{:});

    seeds      = p.Results.Seeds;
    lengths    = sort(p.Results.Lengths);   % ensure ascending order
    outputFile = p.Results.OutputFile;

    % SNR levels for noise corruption (in dB)
    snr_levels = [40, 30, 20, 10];

    % Number of independent realizations
    nRealizations = numel(seeds);

    fprintf('==========================================================\n');
    fprintf('  LyPE Benchmark Dataset Generator\n');
    fprintf('==========================================================\n');
    fprintf('  Seeds:        [%s]\n', num2str(seeds));
    fprintf('  Lengths:      [%s]\n', num2str(lengths));
    fprintf('  SNR levels:   [%s] dB\n', num2str(snr_levels));
    fprintf('  Realizations: %d\n', nRealizations);
    fprintf('  Output file:  %s\n', outputFile);
    fprintf('==========================================================\n\n');

    %% ====================================================================
    %  DATASET A: Van der Pol Oscillator (Periodic / Limit Cycle)
    % =====================================================================
    %
    %  The Van der Pol oscillator is described by:
    %      x'' - mu*(1 - x^2)*x' + x = 0
    %
    %  Rewritten as a first-order system:
    %      dx/dt = y
    %      dy/dt = mu*(1 - x^2)*y - x
    %
    %  With mu = 1, the system exhibits a stable limit cycle.
    %  Theoretical largest Lyapunov exponent: 0 (periodic orbit).
    %
    % =====================================================================
    fprintf('[Dataset A] Van der Pol oscillator (periodic, mu=1)...\n');

    vdp_params.mu           = 1;           % nonlinearity parameter
    vdp_params.dt           = 0.01;        % sampling interval (model time units)
    vdp_params.transient    = 500;         % transient duration to discard (time units)
    vdp_params.ic           = [0.1; 0.0];  % initial conditions [x; x']

    % Maximum number of post-transient points needed
    max_len = max(lengths);

    % Total integration time: transient + enough time for max_len points
    % Each point is dt apart, so we need max_len * dt time units of data
    total_time_needed = vdp_params.transient + max_len * vdp_params.dt;

    % ODE45 options with tight tolerances for accuracy
    ode_opts_vdp = odeset('RelTol', 1e-10, 'AbsTol', 1e-12, ...
                          'MaxStep', vdp_params.dt);

    % Van der Pol ODE right-hand side (defined as nested function below)
    % dx/dt = y
    % dy/dt = mu*(1 - x^2)*y - x
    mu_val = vdp_params.mu;
    vdp_rhs = @(t, state) [state(2); ...
                            mu_val * (1 - state(1)^2) * state(2) - state(1)];

    % Pre-allocate storage
    periodic_clean = cell(nRealizations, 1);

    for iReal = 1:nRealizations
        rng(seeds(iReal), 'twister');
        fprintf('  Realization %d/%d (seed=%d)...', iReal, nRealizations, seeds(iReal));

        % Slightly perturb initial conditions for each realization
        % (all should converge to the same limit cycle)
        if iReal == 1
            ic = vdp_params.ic;
        else
            ic = vdp_params.ic + 0.01 * randn(2, 1);
        end

        % Integration time span with fine output sampling
        tspan = 0 : vdp_params.dt : total_time_needed;

        % Integrate the Van der Pol oscillator
        [t_full, y_full] = ode45(vdp_rhs, tspan, ic, ode_opts_vdp);

        % Determine transient cutoff index
        transient_idx = find(t_full >= vdp_params.transient, 1, 'first');

        % Extract post-transient data (x-component only for univariate analysis)
        y_post = y_full(transient_idx:end, :);

        % Store data at each requested length
        periodic_clean{iReal} = cell(numel(lengths), 1);
        for iLen = 1:numel(lengths)
            N = lengths(iLen);
            if size(y_post, 1) < N
                error('generate_benchmarks:insufficientData', ...
                    'Not enough post-transient VdP data for length %d. Got %d points.', ...
                    N, size(y_post, 1));
            end
            % Store x-component (column 1) as the primary time series
            periodic_clean{iReal}{iLen} = y_post(1:N, 1);
        end
        fprintf(' done.\n');
    end

    % Generate noise-corrupted versions
    fprintf('  Generating noisy versions (SNR = %s dB)...\n', num2str(snr_levels));
    periodic_noisy = cell(numel(snr_levels), 1);
    for iSNR = 1:numel(snr_levels)
        periodic_noisy{iSNR} = cell(nRealizations, 1);
        for iReal = 1:nRealizations
            rng(seeds(iReal) + 100 * iSNR, 'twister');  % reproducible noise
            periodic_noisy{iSNR}{iReal} = cell(numel(lengths), 1);
            for iLen = 1:numel(lengths)
                x_clean = periodic_clean{iReal}{iLen};
                noise_std = std(x_clean) / (10^(snr_levels(iSNR) / 20));
                x_noisy = x_clean + noise_std * randn(size(x_clean));
                periodic_noisy{iSNR}{iReal}{iLen} = x_noisy;
            end
        end
    end

    % Assemble the periodic benchmark struct
    benchmarks.periodic.clean          = periodic_clean;
    benchmarks.periodic.noisy          = periodic_noisy;
    benchmarks.periodic.params         = vdp_params;
    benchmarks.periodic.theoretical_lye = 0;
    benchmarks.periodic.time_unit      = 'model_time_unit';
    benchmarks.periodic.seeds          = seeds;
    benchmarks.periodic.lengths        = lengths;
    benchmarks.periodic.snr_levels     = snr_levels;
    benchmarks.periodic.description    = ['Van der Pol oscillator (mu=1). ' ...
        'Stable limit cycle with theoretical LyE = 0. ' ...
        'Data is the x-component sampled at dt=0.01 model time units.'];

    fprintf('  [Dataset A] Complete.\n\n');

    %% ====================================================================
    %  DATASET B: Lorenz-63 Chaotic System
    % =====================================================================
    %
    %  The Lorenz system (1963):
    %      dx/dt = sigma * (y - x)
    %      dy/dt = x * (rho - z) - y
    %      dz/dt = x * y - beta * z
    %
    %  Standard chaotic parameters: sigma=10, rho=28, beta=8/3
    %  Reference largest Lyapunov exponent: ~0.9056 per model time unit
    %  (Sprott 2003, Table A.1; Viswanath 1998)
    %
    % =====================================================================
    fprintf('[Dataset B] Lorenz-63 chaotic system...\n');

    lorenz_params.sigma     = 10;
    lorenz_params.rho       = 28;
    lorenz_params.beta      = 8/3;
    lorenz_params.dt        = 0.01;        % sampling interval (model time units)
    lorenz_params.transient = 1000;        % transient to discard (time units)
    lorenz_params.ic        = [1; 1; 1];   % initial conditions [x; y; z]

    % Total integration time
    total_time_lorenz = lorenz_params.transient + max_len * lorenz_params.dt;

    % ODE45 options — very tight tolerances for chaotic system accuracy
    ode_opts_lorenz = odeset('RelTol', 1e-12, 'AbsTol', 1e-14, ...
                             'MaxStep', lorenz_params.dt);

    % Lorenz ODE right-hand side
    sigma_val = lorenz_params.sigma;
    rho_val   = lorenz_params.rho;
    beta_val  = lorenz_params.beta;
    lorenz_rhs = @(t, s) [sigma_val * (s(2) - s(1)); ...
                           s(1) * (rho_val - s(3)) - s(2); ...
                           s(1) * s(2) - beta_val * s(3)];

    % Pre-allocate storage: store x, y, z components separately
    lorenz_clean_x = cell(nRealizations, 1);
    lorenz_clean_y = cell(nRealizations, 1);
    lorenz_clean_z = cell(nRealizations, 1);

    for iReal = 1:nRealizations
        rng(seeds(iReal), 'twister');
        fprintf('  Realization %d/%d (seed=%d)...', iReal, nRealizations, seeds(iReal));

        % Perturb initial conditions slightly for different realizations
        if iReal == 1
            ic = lorenz_params.ic;
        else
            ic = lorenz_params.ic + 0.001 * randn(3, 1);
        end

        % Integration time span
        tspan = 0 : lorenz_params.dt : total_time_lorenz;

        % Integrate the Lorenz system
        [t_full, y_full] = ode45(lorenz_rhs, tspan, ic, ode_opts_lorenz);

        % Discard transient
        transient_idx = find(t_full >= lorenz_params.transient, 1, 'first');
        y_post = y_full(transient_idx:end, :);

        % Store each component at each requested length
        lorenz_clean_x{iReal} = cell(numel(lengths), 1);
        lorenz_clean_y{iReal} = cell(numel(lengths), 1);
        lorenz_clean_z{iReal} = cell(numel(lengths), 1);

        for iLen = 1:numel(lengths)
            N = lengths(iLen);
            if size(y_post, 1) < N
                error('generate_benchmarks:insufficientData', ...
                    'Not enough post-transient Lorenz data for length %d. Got %d points.', ...
                    N, size(y_post, 1));
            end
            lorenz_clean_x{iReal}{iLen} = y_post(1:N, 1);  % x-component
            lorenz_clean_y{iReal}{iLen} = y_post(1:N, 2);  % y-component
            lorenz_clean_z{iReal}{iLen} = y_post(1:N, 3);  % z-component
        end
        fprintf(' done.\n');
    end

    % Generate noise-corrupted versions (applied independently to each component)
    fprintf('  Generating noisy versions (SNR = %s dB)...\n', num2str(snr_levels));
    lorenz_noisy_x = cell(numel(snr_levels), 1);
    lorenz_noisy_y = cell(numel(snr_levels), 1);
    lorenz_noisy_z = cell(numel(snr_levels), 1);

    for iSNR = 1:numel(snr_levels)
        lorenz_noisy_x{iSNR} = cell(nRealizations, 1);
        lorenz_noisy_y{iSNR} = cell(nRealizations, 1);
        lorenz_noisy_z{iSNR} = cell(nRealizations, 1);

        for iReal = 1:nRealizations
            lorenz_noisy_x{iSNR}{iReal} = cell(numel(lengths), 1);
            lorenz_noisy_y{iSNR}{iReal} = cell(numel(lengths), 1);
            lorenz_noisy_z{iSNR}{iReal} = cell(numel(lengths), 1);

            for iLen = 1:numel(lengths)
                % Use different RNG offsets for x, y, z to get independent noise
                rng(seeds(iReal) + 200 * iSNR + 1000 * iLen, 'twister');

                x_clean = lorenz_clean_x{iReal}{iLen};
                y_clean = lorenz_clean_y{iReal}{iLen};
                z_clean = lorenz_clean_z{iReal}{iLen};

                noise_std_x = std(x_clean) / (10^(snr_levels(iSNR) / 20));
                noise_std_y = std(y_clean) / (10^(snr_levels(iSNR) / 20));
                noise_std_z = std(z_clean) / (10^(snr_levels(iSNR) / 20));

                lorenz_noisy_x{iSNR}{iReal}{iLen} = x_clean + noise_std_x * randn(size(x_clean));
                lorenz_noisy_y{iSNR}{iReal}{iLen} = y_clean + noise_std_y * randn(size(y_clean));
                lorenz_noisy_z{iSNR}{iReal}{iLen} = z_clean + noise_std_z * randn(size(z_clean));
            end
        end
    end

    % Assemble the Lorenz benchmark struct
    benchmarks.lorenz.clean.x          = lorenz_clean_x;
    benchmarks.lorenz.clean.y          = lorenz_clean_y;
    benchmarks.lorenz.clean.z          = lorenz_clean_z;
    benchmarks.lorenz.noisy.x          = lorenz_noisy_x;
    benchmarks.lorenz.noisy.y          = lorenz_noisy_y;
    benchmarks.lorenz.noisy.z          = lorenz_noisy_z;
    benchmarks.lorenz.params           = lorenz_params;
    benchmarks.lorenz.theoretical_lye  = 0.9056;
    benchmarks.lorenz.time_unit        = 'model_time_unit';
    benchmarks.lorenz.seeds            = seeds;
    benchmarks.lorenz.lengths          = lengths;
    benchmarks.lorenz.snr_levels       = snr_levels;
    benchmarks.lorenz.description      = ['Lorenz-63 system (sigma=10, rho=28, beta=8/3). ' ...
        'Chaotic attractor with reference LyE ~ 0.9056 per model time unit ' ...
        '(Sprott 2003, Viswanath 1998). Components x, y, z stored separately ' ...
        'for univariate and multivariate analysis. Sampled at dt=0.01.'];

    fprintf('  [Dataset B] Complete.\n\n');

    %% ====================================================================
    %  DATASET C: Logistic Map (Fully Chaotic Regime)
    % =====================================================================
    %
    %  The logistic map:
    %      x(n+1) = r * x(n) * (1 - x(n))
    %
    %  With r = 4, the map is fully chaotic on [0, 1].
    %  Theoretical largest Lyapunov exponent: ln(2) ≈ 0.6931 per iteration.
    %
    %  This is an exact result. The invariant measure for r=4 is the
    %  arcsine distribution, and the Lyapunov exponent can be computed
    %  analytically via the change of variables x = sin^2(pi*theta).
    %
    % =====================================================================
    fprintf('[Dataset C] Logistic map (r=4, fully chaotic)...\n');

    logistic_params.r         = 4;
    logistic_params.transient = 1000;      % transient iterates to discard
    logistic_params.x0        = 0.1;       % initial condition (avoid 0, 0.25, 0.5, 0.75, 1)

    % Total iterates needed: transient + max data length
    total_iterates = logistic_params.transient + max_len;

    % Pre-allocate
    logistic_clean = cell(nRealizations, 1);

    for iReal = 1:nRealizations
        rng(seeds(iReal), 'twister');
        fprintf('  Realization %d/%d (seed=%d)...', iReal, nRealizations, seeds(iReal));

        % Perturb initial condition slightly for each realization
        if iReal == 1
            x0 = logistic_params.x0;
        else
            x0 = logistic_params.x0 + 0.001 * randn();
            % Clamp to valid range (0, 1), avoiding fixed points
            x0 = max(0.001, min(0.999, x0));
        end

        % Iterate the logistic map
        x_map = zeros(total_iterates + 1, 1);
        x_map(1) = x0;
        for n = 1:total_iterates
            x_map(n+1) = logistic_params.r * x_map(n) * (1 - x_map(n));
        end

        % Discard transient
        x_post = x_map(logistic_params.transient + 1 : end);

        % Store at each requested length
        logistic_clean{iReal} = cell(numel(lengths), 1);
        for iLen = 1:numel(lengths)
            N = lengths(iLen);
            if numel(x_post) < N
                error('generate_benchmarks:insufficientData', ...
                    'Not enough post-transient logistic data for length %d.', N);
            end
            logistic_clean{iReal}{iLen} = x_post(1:N);
        end
        fprintf(' done.\n');
    end

    % Generate noise-corrupted versions
    fprintf('  Generating noisy versions (SNR = %s dB)...\n', num2str(snr_levels));
    logistic_noisy = cell(numel(snr_levels), 1);
    for iSNR = 1:numel(snr_levels)
        logistic_noisy{iSNR} = cell(nRealizations, 1);
        for iReal = 1:nRealizations
            rng(seeds(iReal) + 300 * iSNR, 'twister');  % reproducible noise
            logistic_noisy{iSNR}{iReal} = cell(numel(lengths), 1);
            for iLen = 1:numel(lengths)
                x_clean = logistic_clean{iReal}{iLen};
                noise_std = std(x_clean) / (10^(snr_levels(iSNR) / 20));
                x_noisy = x_clean + noise_std * randn(size(x_clean));
                logistic_noisy{iSNR}{iReal}{iLen} = x_noisy;
            end
        end
    end

    % Assemble the logistic benchmark struct
    benchmarks.logistic.clean          = logistic_clean;
    benchmarks.logistic.noisy          = logistic_noisy;
    benchmarks.logistic.params         = logistic_params;
    benchmarks.logistic.theoretical_lye = log(2);  % ln(2) ≈ 0.6931
    benchmarks.logistic.time_unit      = 'iteration';
    benchmarks.logistic.seeds          = seeds;
    benchmarks.logistic.lengths        = lengths;
    benchmarks.logistic.snr_levels     = snr_levels;
    benchmarks.logistic.description    = ['Logistic map x(n+1) = 4*x(n)*(1-x(n)). ' ...
        'Fully chaotic regime with theoretical LyE = ln(2) per iteration (exact). ' ...
        'Initial condition x(1) = 0.1, 1000 transient iterates discarded.'];

    fprintf('  [Dataset C] Complete.\n\n');

    %% ====================================================================
    %  DATASET D: Linear Contractive Map (Supplementary)
    % =====================================================================
    %
    %  The linear contractive map:
    %      x(n+1) = a * x(n)     with a = 0.95
    %
    %  This is a trivially simple system with:
    %      x(n) = a^n * x(1) = 0.95^n
    %
    %  Theoretical Lyapunov exponent: ln(a) = ln(0.95) ≈ -0.0513 per iteration.
    %
    %  NOTE: This system converges to a fixed point (x=0) and does NOT
    %  produce a strange attractor. It is NOT suitable for attractor-based
    %  methods like Rosenstein's algorithm. However, it is useful for:
    %    - Verifying distance/log/slope computations
    %    - Testing that algorithms correctly handle contracting dynamics
    %    - Debugging divergence curve calculations
    %
    % =====================================================================
    fprintf('[Dataset D] Linear contractive map (a=0.95, supplementary)...\n');

    contractive_params.a   = 0.95;
    contractive_params.x0  = 1.0;
    contractive_params.len = 4000;  % fixed length for this simple system

    % Pre-allocate
    contractive_clean = cell(nRealizations, 1);

    for iReal = 1:nRealizations
        rng(seeds(iReal), 'twister');
        fprintf('  Realization %d/%d (seed=%d)...', iReal, nRealizations, seeds(iReal));

        % For the contractive map, the deterministic part is identical for
        % all realizations (no dependence on random IC perturbation needed),
        % but we still generate separate realizations for consistency.
        % Add tiny perturbation to x0 for non-trivial realization differences.
        if iReal == 1
            x0 = contractive_params.x0;
        else
            x0 = contractive_params.x0 + 1e-6 * randn();
        end

        N = contractive_params.len;
        x_map = zeros(N, 1);
        x_map(1) = x0;
        for n = 1:N-1
            x_map(n+1) = contractive_params.a * x_map(n);
        end

        % Store as single length (no multiple lengths for this dataset)
        contractive_clean{iReal} = x_map;
        fprintf(' done.\n');
    end

    % Generate noise-corrupted versions for the contractive map
    fprintf('  Generating noisy versions (SNR = %s dB)...\n', num2str(snr_levels));
    contractive_noisy = cell(numel(snr_levels), 1);
    for iSNR = 1:numel(snr_levels)
        contractive_noisy{iSNR} = cell(nRealizations, 1);
        for iReal = 1:nRealizations
            rng(seeds(iReal) + 400 * iSNR, 'twister');
            x_clean = contractive_clean{iReal};
            noise_std = std(x_clean) / (10^(snr_levels(iSNR) / 20));
            x_noisy = x_clean + noise_std * randn(size(x_clean));
            contractive_noisy{iSNR}{iReal} = x_noisy;
        end
    end

    % Assemble the contractive benchmark struct
    benchmarks.contractive.clean          = contractive_clean;
    benchmarks.contractive.noisy          = contractive_noisy;
    benchmarks.contractive.params         = contractive_params;
    benchmarks.contractive.theoretical_lye = log(0.95);  % ln(0.95) ≈ -0.0513
    benchmarks.contractive.time_unit      = 'iteration';
    benchmarks.contractive.seeds          = seeds;
    benchmarks.contractive.lengths        = contractive_params.len;
    benchmarks.contractive.snr_levels     = snr_levels;
    benchmarks.contractive.description    = ['Linear contractive map x(n+1) = 0.95*x(n). ' ...
        'Theoretical LyE = ln(0.95) per iteration (exact). ' ...
        'NOT suitable for attractor-based Rosenstein testing — converges to ' ...
        'fixed point. Useful for verifying distance/log/slope computations.'];

    fprintf('  [Dataset D] Complete.\n\n');

    %% ====================================================================
    %  ADD METADATA AND SAVE
    % =====================================================================
    fprintf('Assembling metadata and saving...\n');

    % Top-level metadata
    benchmarks.metadata.generated_by   = 'generate_benchmarks.m';
    benchmarks.metadata.generated_date = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    benchmarks.metadata.matlab_version = version;
    benchmarks.metadata.seeds          = seeds;
    benchmarks.metadata.lengths        = lengths;
    benchmarks.metadata.snr_levels     = snr_levels;
    benchmarks.metadata.nRealizations  = nRealizations;
    benchmarks.metadata.noise_model    = 'Additive white Gaussian noise: x_noisy = x + noise_std * randn(size(x)), noise_std = std(x) / 10^(SNR_dB/20)';

    % Summary table of theoretical/reference LyE values
    benchmarks.metadata.summary = { ...
        'Dataset',      'System',             'LyE (theory/ref)',    'Time Unit'; ...
        'periodic',     'Van der Pol (mu=1)',  '0',                  'model_time_unit'; ...
        'lorenz',       'Lorenz-63',           '~0.9056',            'model_time_unit'; ...
        'logistic',     'Logistic (r=4)',      'ln(2) ~ 0.6931',     'iteration'; ...
        'contractive',  'Linear (a=0.95)',     'ln(0.95) ~ -0.0513', 'iteration'};

    % Save to MAT file
    save(outputFile, 'benchmarks', '-v7.3');
    fprintf('\nBenchmark data saved to: %s\n', outputFile);

    % Print summary
    fprintf('\n==========================================================\n');
    fprintf('  BENCHMARK GENERATION SUMMARY\n');
    fprintf('==========================================================\n');
    fprintf('  %-15s  %-12s  %-18s  %s\n', 'Dataset', 'LyE', 'Time Unit', 'Lengths');
    fprintf('  %-15s  %-12s  %-18s  %s\n', '-------', '---', '---------', '-------');
    fprintf('  %-15s  %-12.4f  %-18s  %s\n', 'periodic',     benchmarks.periodic.theoretical_lye,     benchmarks.periodic.time_unit,     num2str(lengths));
    fprintf('  %-15s  %-12.4f  %-18s  %s\n', 'lorenz',       benchmarks.lorenz.theoretical_lye,       benchmarks.lorenz.time_unit,       num2str(lengths));
    fprintf('  %-15s  %-12.4f  %-18s  %s\n', 'logistic',     benchmarks.logistic.theoretical_lye,     benchmarks.logistic.time_unit,     num2str(lengths));
    fprintf('  %-15s  %-12.4f  %-18s  %s\n', 'contractive',  benchmarks.contractive.theoretical_lye,  benchmarks.contractive.time_unit,  num2str(contractive_params.len));
    fprintf('==========================================================\n');
    fprintf('  Each dataset: %d realizations x %d SNR levels + clean\n', nRealizations, numel(snr_levels));
    fprintf('  Total file size: %.1f MB\n', dir(outputFile).bytes / 1e6);
    fprintf('==========================================================\n');
    fprintf('Done.\n');

end
