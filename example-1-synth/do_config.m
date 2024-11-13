% Synthetic Field Trip data analysis - Configuration
% Written by Christopher Thomas.


%
% Switches that affect static configuration.

% Number of phase bins.
phase_bin_density = 'small';
%phase_bin_density = 'medium';
%phase_bin_density = 'large';

% Number of frequency bins.
freq_bin_density = 'small';
%freq_bin_density = 'medium';
%freq_bin_density = 'large';

% Cross-correlation window time step granularity.
%xcorr_time_step_ms = 50;
%xcorr_time_step_ms = 200;
xcorr_time_step_ms = 500;

% Cross-correlation time lag granularity.
% Set this to 0 for a step size of one sample.
%xcorr_delay_step_ms = 0;
%xcorr_delay_step_ms = 1;
xcorr_delay_step_ms = 2;

% Extrapolation for MI and TE calculations.
want_mi_extrap = false;
want_te_extrap = false;

% Reduced bin counts for more accurate entropy estimates.
want_fewer_bins = true;

want_replicates = false;
want_fewer_replicates = true;



%
% Configuration that rarely changes.

do_config_static;



%
% Configuration that changes frequently.


% Dataset to use.

%dataset = dataset_smoke;
%dataset = dataset_short;
dataset = dataset_ab;
%dataset = dataset_loop;


% Filtering to use.
% band_ [alpha/beta/gamma] _ [full/narrow] (or band_none to not filter).
%bandpass_band = band_none;
bandpass_band = band_alpha_full;


% Multithreading. Requires Parallel Computing Toolbox.

want_parallel = true;


% Types of plot to generate.

want_wave_plots = false;
want_timelock_plots = false;
want_spect_plots = false;

want_phase_stats_report = false;

want_do_xc = true;
want_do_pc = true;
want_do_mi = false;
want_do_te = false;

% Merge trials instead of averaging across them.
want_span_trials = true;

want_info_plots_lagtime = true;
want_info_plots_line = true;
want_info_plots_search = false;
want_info_plots_stats_box = true;
want_info_plots_stats_grid = true;
want_info_plots_stats_report = true;
want_info_plots_network = true;
want_info_plots_network_anim = false;

want_info_plots_byphase = false;
want_info_plots_byphase_more = false;
want_info_plots_byfreq = false;

want_htplots_coherence = false;
want_htplots_powercorrel = false;
want_htplots_powernongauss = false;


% Waveform plot config.

wave_plots_wanted = {};
wave_plots_wanted = [ wave_plots_wanted { 'oneplot' } ];
%wave_plots_wanted = [ wave_plots_wanted { 'pertrial' } ];
%wave_plots_wanted = [ wave_plots_wanted { 'perchannel' } ];

wave_zoom_ranges = { [] };
wave_zoom_labels = { 'full' };


% Timelock plot config.

timelock_plots_wanted = {};
timelock_plots_wanted = [ timelock_plots_wanted { 'oneplot' } ];
%timelock_plots_wanted = [ timelock_plots_wanted { 'perchannel' } ];
%timelock_plots_wanted = [ timelock_plots_wanted { 'stripchart' } ];

timelock_zoom_ranges_ms = { [] };
timelock_zoom_labels = { 'full' };


% Spectrogram plot config.

spect_zoom_ranges = { [], [ -1 2 ] };
spect_zoom_labels = { 'all', 'close' };

spect_plot_freqrange = [];
spect_loglin = 'linear';


% Cross-correlation plot config.

% Nothing yet.


%
% This is the end of the file.
