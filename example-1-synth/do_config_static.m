% Synthetic Field Trip data analysis - Configuration that rarely changes.
% Written by Christopher Thomas.


%
% Folders.

plotdir = 'plots';

datadir = [ '..' filesep 'datasets' filesep '20231214-synth-10tr' ];



%
% Signal processing configuration.

% NOTE - The lowest frequency should fit in the window!
spectrogram_win_secs = 0.5;
spectrogram_step_secs = 0.02;
spectrogram_freq_range = [ 2 100 ];


% Different narrow-band options.

% "don't filter" option.
band_none = [];

% Full standard bands.
band_alpha_full = [ 6 12 ];
band_beta_full = [ 12 24 ];
band_gamma_full = [ 25 50 ];

% Trying to extract the 8 Hz peak and its harmonics.
band_alpha_narrow = [ 7 9 ];
band_beta_narrow = [ 14 18 ];
band_gamma_narrow = [ 28 36 ];


%
% FIXME - The "xcorr" stuff is for all information measures, not just
% cross-correlation. Renaming it would be a massive pain.

% This is 'detrend', 'zeromean', or 'none'.
xcorr_detrend = 'detrend';

% Specify whether we want r or r^2 for Pearson's.
% MI for Gaussian RVs is -(1/2) log_2( 1 - r^2 ).
% r looks almost identical to XC, r^2 looks qualitatively like MI.
pcorr_want_r2 = true;

% NOTE - The system's response time is on the order of 1 second, so we need
% a much longer comparison time!
xcorr_want_long = true;

% This is the "win_params" structure the helper function wants, with extra
% fields. NOTE - Formerly "mua_params".
xcorr_params = struct();

% Lag range. Nominal delays are < 10 ms. Osc period for alpha is 120.
xcorr_params.delay_range_ms = [ -100 100 ];
%xcorr_params.delay_range_ms = [ -50 50 ];
xcorr_params.delay_step_ms = xcorr_delay_step_ms;

% This is 'unbiased' to normalize by sample count, and 'coeff' to make
% self-correlation equal to 1.
xcorr_norm_method = 'coeff';

if xcorr_want_long
  % NOTE - This will have a large memory footprint if we have 128ch!

  xcorr_params.time_window_ms = 500;

  % Long dataset ran from -10 to +15.
%  xcorr_params.timelist_ms = [ -5000:xcorr_time_step_ms:10000 ];

  % Short dataset ran from -5 to +10.
  xcorr_params.timelist_ms = [ -4000:xcorr_time_step_ms:8000 ];
else
  xcorr_params.time_window_ms = 300;

  % NOTE - We don't need "time_before_ms", since xcorr only cares about
  % timelist_ms.

  % For real data, we normally have a gap around stimulation, but for the
  % synthetic data we can just have a uniform span.

  %xcorr_params.timelist_ms = ...
    [ (-500:xcorr_time_step_ms:-200) (200:xcorr_time_step_ms:1000) ];
  xcorr_params.timelist_ms = [ -500:xcorr_time_step_ms:1000 ];
end

% FIXME - We want smaller time windows for phase-binned XC, since alignment
% is transient.
xcorr_params_phase = xcorr_params;
xcorr_params_phase.time_window_ms = 300;


if want_fewer_bins
mutual_bins = 4;
transfer_bins = 2;
else
mutual_bins = 8;
%transfer_bins = 8;
transfer_bins = 4;
end


pcorr_replicates = 0;
mutual_replicates = 0;
transfer_replicates = 0;

if want_replicates
  if want_fewer_replicates
    pcorr_replicates = 10;
    mutual_replicates = 3;
    transfer_replicates = 3;
  else
    pcorr_replicates = 30;
    mutual_replicates = 10;
    transfer_replicates = 10;
  end
end


% Set this to struct() for default extrapolation, 'none' for no extrapolation.

mutual_extrap = 'none';
transfer_extrap = 'none';

if want_mi_extrap
  mutual_extrap = struct();
end
if want_te_extrap
  transfer_extrap = struct();
end


% Various other xcorr analysis config not stored in the parameter structure.

if xcorr_want_long
  xcorr_after_start_ms = 1000;
  xcorr_before_start_ms = 0;
else
  xcorr_after_start_ms = 200;
  xcorr_before_start_ms = 0;
end

if strcmp('large', phase_bin_density)
  % High-density analysis. 12 degree intervals.
  xcorr_phase_bins = 30;
  xcorr_phase_bin_width_deg = 24;
elseif strcmp('medium', phase_bin_density)
  % Medium-density analysis. 45 degree intervals.
  xcorr_phase_bins = 8;
  xcorr_phase_bin_width_deg = 60;
else
  % Fast analysis for testing. 120 degree intervals.
  xcorr_phase_bins = 3;
  xcorr_phase_bin_width_deg = 120;
end


% Frequency bands, for cross-correlation vs frequency.

if strcmp('large', freq_bin_density)
  % Octave bands with 75% overlap.
  xcorr_vsfreq_bands = ...
    { [2.8 5.6], [3.4 6.8], [4.0 8.0], [4.8 9.6], ...
      [5.6 11], [6.8 14], [8 16], [9.6 19], ...
      [11 22], [14 28], [16 32], [19 38], ...
      [22 44], [28 56], [32 64], [38 76], ...
      [44 88], [56 112] };
elseif strcmp('medium', freq_bin_density)
  % Octave bands with 50% overlap.
  xcorr_vsfreq_bands = ...
    { [3 6], [4 8], [6 12], [8 16], [12 24], [16 32], ...
      [24 48], [32 64], [48 96] };
else
  % A small number of bands for testing.
  xcorr_vsfreq_bands = { [6 12], [8 16], [12 24] };
end


% Time smoothing window for collapsing time for phase-binned plots.
%xcorr_phase_time_smooth_ms = NaN;
%xcorr_phase_time_smooth_ms = 200;
xcorr_phase_time_smooth_ms = 500;
%xcorr_phase_time_smooth_ms = 1000;

% Minimum XC amplitude for contributing to phase-binned averaging.
% This was intended as a proxy for phase locking, but doesn't seem to
% make much difference, and we can directly filter on PLV.
xcorr_phase_min_magnitude = -1;
%xcorr_phase_min_magnitude = 0.05;

% Min phase-lock value for phase binning.
% PLV is 1 - circvar.

% FIXME - These values were computed incorrectly. Actual is 0.8-0.95.
% 500 ms windows gives average/median PLV of about 0.05 +/- 0.03 or so.
% 300 ms windows gives 0.11 +/- 0.03 or so.

% If an absolute threshold is defined, it'll take precedence over a
% percentile threshold. Use NaN to fall back to percentile.
% NaN for both means "don't filter on PLV".

% Absolute PLV threshold.
%xcorr_phase_min_plv = NaN;
xcorr_phase_min_plv = 0.5;
%xcorr_phase_min_plv = 0.7;

% Percentile PLV threshold.
%xcorr_phase_plv_prctile = NaN;
%xcorr_phase_plv_prctile = 50;
xcorr_phase_plv_prctile = 10;


% Time smoothing window for collapsing lag/time.
% NaN to not smooth.
%xcorr_collapse_smooth_ms = 1000;
xcorr_collapse_smooth_ms = 3000;


% Time smoothing window for peak detection.
% NaN to not smooth.

% NOTE - Changing this to be a multiple of midband period.
% Original 6-12 Hz smoothing was 1/2/5 sec (8/16/40 wavelengths).
xcorr_peak_smooth_wavecount = 20;
%xcorr_peak_smooth_wavecount = 40;

if true
  % Method used for identifying nominal xcorr peak.
  xcorr_peak_method = 'largest';
  % Lag search range for the largest peak; [] for full range.
%  xcorr_peak_location = [];
  xcorr_peak_location = [ -30 30 ];
else
  % Method used for identifying nominal xcorr peak.
  xcorr_peak_method = 'nearest';
  % Starting lag value for finding the nearest peak.
  xcorr_peak_location = 0;
end

% Amplitude threshold needed for time-series peaks to be valid.
% Set this NaN to disable filtering.
xcorr_time_amp_threshold = NaN;
%xcorr_time_amp_threshold = 0.04;
%xcorr_time_amp_threshold = 0.06;


% Tuning parameters for black magic cross-correlation statistics.

% Filtering to decide what lag range to search and what amplitudes to accept.
xcorr_peak_halo_thresh = 0;
xcorr_peak_mag_accept = [ 0.5 inf ];

xcorr_stats_config = struct();

% Reject any pair cross-correlation magnitude that varies too much.
%xcorr_stats_config.max_amp_dev = 0.01;
%xcorr_stats_config.max_amp_dev = 0.015;
xcorr_stats_config.max_amp_dev = 0.02;
%xcorr_stats_config.max_amp_dev = inf;

% Reject any pair lag that varies too much.
%xcorr_stats_config.max_lag_dev = 2;
xcorr_stats_config.max_lag_dev = 3;

% Reject any pair lag for pairs that had poor correlation.
% FIXME - This has to be a relative measure for MI and TE!
%xcorr_stats_config.min_amp_for_lag = 0.05;
xcorr_stats_config.min_amp_for_lag = -1;
% NOTE - For XC, setting this to 0.2 loses a lot. Assume max XC _is_ 1.
%xcorr_stats_config.min_amp_rel_for_lag = 0.2;
xcorr_stats_config.min_amp_rel_for_lag = 0.05;

% Statistical significance thresholds for changes (sigma values).
xcorr_stats_config.significant_amp_sigma = 2;
xcorr_stats_config.significant_lag_sigma = 2;

% Whether to squash all statistics of a channel pair if any of them are bad.
xcorr_squash_if_any_nan = false;


% Tuning parameters for plotting animated cross-correlation network graphs.

xcorr_network_config = struct();

xcorr_network_config.amp_range = [ 0.03 0.20 ];
xcorr_network_config.lag_range_ms = [ 5 20 ];

xcorr_network_config.amp_rolloff = { 'nan', 'sigmoid' };
xcorr_network_config.lag_rolloff = { 'nan', 'clamp' };

xcorr_network_config.layout = 'circle';
%xcorr_network_config.layout = 'adaptive';


% Per-case tuning for network graphs.

% NOTE - 'auto' normalizes each case separately, so it's not very useful.
xcorr_network_config.amp_auto_prctiles = [ 50 90 ];

xcorr_network_amp_tuning = struct();
xcorr_network_amp_tuning.xc = [ 0.03 0.20 ];

% MI with 8 bins.
%xcorr_network_amp_tuning.mi = [ 0.010 0.030 ];

% MI with 4 bins.
% Auto-tuned was 50/90% 0.007 / 0.012 for most.
xcorr_network_amp_tuning.mi = [ 0.010 0.025 ];

% MI auto.
%xcorr_network_amp_tuning.mi = 'auto';

% TE with 8 bins.
% Auto-tuned 10/50/90 was 0.041 / 0.042 / 0.044 (up 0.048).
%xcorr_network_amp_tuning.te = [ 0.042 0.050 ];
%xcorr_network_amp_tuning.te = [ 0.043 0.050 ];

% TE with 4 bins.
% Auto was 0.008-0.009 median, 0.010 90th for most, 0.015-0.016 for loop.
%xcorr_network_amp_tuning.te = [ 0.0085 0.0150 ];

% TE with 2 bins.
% Auto was about 0.0016 median, 0.0025 90th for most, 0.0011-0.0014 loopdown.
xcorr_network_amp_tuning.te = [ 0.0020 0.0040 ];

% TE auto.
%xcorr_network_amp_tuning.te = 'auto';



%
% Plotting configuration.

% Yellow/cyan colour map, for amplitude.
colmapamp = nlPlot_getColorMapHotCold( ...
  [ 0.6 0.9 1.0 ], [ 1.0 0.9 0.3 ], 1.0 );

% Green/magenta colour map, for time lag.
colmaplag = nlPlot_getColorMapHotCold( ...
  [ 0.9 0.3 0.9 ], [ 0.3 0.9 0.3 ], 1.0 );


% These plots aren't particularly useful.
xcorr_plot_want_timemean = false;
xcorr_plot_want_lagdev = false;

% These should be plot-, filename-, and fieldname-safe.
xcorr_plot_list_vslag = {};
%xcorr_plot_list_vslag = [ xcorr_plot_list_vslag {'all'} ];
%xcorr_plot_list_vslag = [ xcorr_plot_list_vslag {'before'} ];
xcorr_plot_list_vslag = [ xcorr_plot_list_vslag {'after'} ];

xcorr_plot_devcursors = [ 0 ];
xcorr_plot_meancursors = [ 0, -0.05, 0.05 ];
xcorr_plot_ampcursors = [ 0, -0.05, 0.05 ];
xcorr_plot_lagcursors = [ 0 ];

% Whether to render NaN as empty (NaN) or as some data value.
%xcorr_plot_phasenan = NaN;
xcorr_plot_phasenan = 0;

% FIXME - Stuff the animation frame rate into the network config structure.
% If we don't supply one, a default rate is used.
xcorr_anim_speedup = 1.0;
xcorr_network_config.anim_framerate = ...
  round(xcorr_anim_speedup * 1000 / xcorr_time_step_ms);



%
% Dataset configurations.

% The synthetic data is called "MUA" but actually looks like LFP, since it's
% just average firing rates rather than the sum of individual spikes.

dataset_short = struct( ...
  'settitle', 'Small Dataset (10 trials)', ...
  'setlabel', 'small', ...
  'paths', {{ ...
    [ datadir filesep 'mua-loopbase-baseline.mat' ], ...
    [ datadir filesep 'mua-loopbase-strong.mat' ], ...
    [ datadir filesep 'mua-loopbase-weak.mat' ], ...
    [ datadir filesep 'mua-loopbase-routeA.mat' ], ...
    [ datadir filesep 'mua-loopbase-routeB.mat' ], ...
    [ datadir filesep 'mua-loopbase-loopup.mat' ], ...
    [ datadir filesep 'mua-loopbase-loopdown.mat' ] }}, ...
  'labels', {{ 'baseline', 'strong', 'weak', ...
    'routeA', 'routeB', 'loopup', 'loopdown' }}, ...
  'titles', {{ 'Baseline', 'Stronger MUA', 'Weaker MUA', ...
    'Routing A', 'Routing B', 'Cycle 1 to 4', 'Cycle 4 to 1' }} );

% The smoke test is a minimum subset of the small dataset.
dataset_smoke = struct( ...
  'settitle', 'Smoke Test (10 trials)', ...
  'setlabel', 'smoke', ...
  'paths', {{ ...
    [ datadir filesep 'mua-loopbase-strong.mat' ] }}, ...
  'labels', {{ 'strong' }}, ...
  'titles', {{ 'Stronger MUA' }} );

% Subsets of the short dataset with less runtime.
dataset_ab = struct( ...
  'settitle', 'A-B Cases (10 trials)', ...
  'setlabel', 'ab', ...
  'paths', {{ ...
    [ datadir filesep 'mua-loopbase-routeA.mat' ], ...
    [ datadir filesep 'mua-loopbase-routeB.mat' ] }}, ...
  'labels', {{ 'routeA', 'routeB' }}, ...
  'titles', {{ 'Routing A', 'Routing B' }} );

dataset_loop = struct( ...
  'settitle', 'Loop Cases (10 trials)', ...
  'setlabel', 'loop', ...
  'paths', {{ ...
    [ datadir filesep 'mua-loopbase-loopup.mat' ], ...
    [ datadir filesep 'mua-loopbase-loopdown.mat' ] }}, ...
  'labels', {{ 'loopup', 'loopdown' }}, ...
  'titles', {{ 'Cycle 1 to 4', 'Cycle 4 to 1' }} );



%
% This is the end of the file.
