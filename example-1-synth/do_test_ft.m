% Synthetic Field Trip data analysis demo script.
% Written by Christopher Thomas.
%
% This is intended to be a short script for analyzing synthetic Field Trip
% data to pick out MUA responses and routing states.
%
% This is used for smoke-testing synthetic data and as an example script for
% doing more thorough analyses.


%
% Initialization.


do_setup_stuff;
do_config;

if want_parallel
  parpool;
end



%
% Build the information analysis cases and function handles.


% Common.

if want_span_trials
  infoflags = { 'spantrials' };
else
  infoflags = { 'avgtrials' };
end

if want_parallel
  infoflags = [ infoflags, { 'parallel' } ];
end


% Cross-correlation.

[ xcpreproc xcparams ] = ...
  eiCalc_getParamsXCorr( xcorr_detrend, xcorr_norm_method );

info_handle_xc = ...
  @( ftdata_dest, ftdata_src, win_params ) ...
  eiCalc_doTimeAndLagAnalysis( ...
    ftdata_dest, ftdata_src, win_params, infoflags, ...
    xcpreproc, @eiCalc_helper_analyzeXCorr, xcparams, ...
    {}, @eiCalc_helper_filterNone, struct() );

info_handle_xc_phase = ...
  @( ftdata_dest, ftdata_src, win_params, phase_params ) ...
  eiCalc_doTimeAndLagAnalysis( ...
    ftdata_dest, ftdata_src, win_params, infoflags, ...
    xcpreproc, @eiCalc_helper_analyzeXCorr, xcparams, ...
    { 'detrend', 'angle', }, @eiCalc_helper_filterPhase, phase_params );


% Pearson's correlation.

pcparams = eiCalc_getParamsPCorr( ...
  pcorr_replicates, pcorr_want_r2, infoflags );

info_handle_pc = ...
  @( ftdata_dest, ftdata_src, win_params ) ...
  eiCalc_doTimeAndLagAnalysis( ...
    ftdata_dest, ftdata_src, win_params, infoflags, ...
    {}, @eiCalc_helper_analyzePCorr, pcparams, ...
    {}, @eiCalc_helper_filterNone, struct() );

info_handle_pc_phase = ...
  @( ftdata_dest, ftdata_src, win_params, phase_params ) ...
  eiCalc_doTimeAndLagAnalysis( ...
    ftdata_dest, ftdata_src, win_params, infoflags, ...
    {}, @eiCalc_helper_analyzePCorr, pcparams, ...
    { 'detrend', 'angle', }, @eiCalc_helper_filterPhase, phase_params );


% Mutual information.

miparams = eiCalc_getParamsMutualTE( ...
  mutual_bins, mutual_bins, mutual_replicates, infoflags, mutual_extrap );

info_handle_mi = ...
  @( ftdata_dest, ftdata_src, win_params ) ...
  eiCalc_doTimeAndLagAnalysis( ...
    ftdata_dest, ftdata_src, win_params, infoflags, ...
    {}, @eiCalc_helper_analyzeMutual, miparams, ...
    {}, @eiCalc_helper_filterNone, struct() );

info_handle_mi_phase = ...
  @( ftdata_dest, ftdata_src, win_params, phase_params ) ...
  eiCalc_doTimeAndLagAnalysis( ...
    ftdata_dest, ftdata_src, win_params, infoflags, ...
    {}, @eiCalc_helper_analyzeMutual, miparams, ...
    { 'detrend', 'angle', }, @eiCalc_helper_filterPhase, phase_params );


% Transfer entropy.

teparams = eiCalc_getParamsMutualTE( ...
  transfer_bins, transfer_bins, transfer_replicates, ...
  infoflags, transfer_extrap );

info_handle_te = ...
  @( ftdata_dest, ftdata_src, win_params ) ...
  eiCalc_doTimeAndLagAnalysis( ...
    ftdata_dest, ftdata_src, win_params, infoflags, ...
    {}, @eiCalc_helper_analyzeTransfer, teparams, ...
    {}, @eiCalc_helper_filterNone, struct() );

info_handle_te_phase = ...
  @( ftdata_dest, ftdata_src, win_params, phase_params ) ...
  eiCalc_doTimeAndLagAnalysis( ...
    ftdata_dest, ftdata_src, win_params, infoflags, ...
    {}, @eiCalc_helper_analyzeTransfer, teparams, ...
    { 'detrend', 'angle', }, @eiCalc_helper_filterPhase, phase_params );


% Nx6 cell array holding tuples of:
% { file label, long name, abbreviation, fieldname, handle, handle_phase }

infocases = {};

if want_do_xc
  infocases = [ infocases ; ...
    { 'xc', 'Cross-Correlation', 'XC', 'xcorr', ...
      info_handle_xc, info_handle_xc_phase } ];
end

if want_do_pc
  infocases = [ infocases ; ...
    { 'pc', 'Pearson''s Correlation', 'PC', 'pcorr', ...
      info_handle_pc, info_handle_pc_phase } ];
end

if want_do_mi
  infocases = [ infocases ; ...
    { 'mi', 'Mutual Information', 'MI', 'mutual', ...
      info_handle_mi, info_handle_mi_phase } ];
end

if want_do_te
  infocases = [ infocases ; ...
    { 'te', 'Transfer Entropy', 'TE', 'transfer', ...
      info_handle_te, info_handle_te_phase } ];
end


infocount = size(infocases,1);



%
% Walk through the dataset's cases, generating plots.


% Make globally aggregated data, too.
% FIXME - Add field names for each info type.
globalphasestatreport = '';
globalxcorrstatreport = '';


disp([ '== Processing "' dataset.settitle '".' ]);

for didx = 1:length(dataset.paths)

  thisfile = dataset.paths{didx};
  thislabel = dataset.labels{didx};
  thistitle = dataset.titles{didx};

  disp([ '-- Loading "' thisfile '".' ]);

  % This loads the following variables:
  % ftdata_mua, modelparams, intcouplings, mixing_matrix, mixing_delays_ms,
  % origheader, origconfigtrl
  % NOTE - We're saving them in a structure, to avoid namespace issues.

  filedata = load(thisfile);

  % Save the raw FT data under a name that we aren't going to reuse for
  % anything else.
  ftdata_wb = filedata.ftdata_mua;

  % Extract metadata.

  samprate_orig = filedata.origheader.Fs;
  origconfigtrl = filedata.origconfigtrl;

  samprate_wb = ftdata_wb.fsample;
  chancount = length(ftdata_wb.label);
  [ chanlabels chantitles ] = euUtil_makeSafeStringArray( ftdata_wb.label );


  %
  % Waveform plots.

  if want_wave_plots
    euPlot_plotFTTrials( ftdata_wb, samprate_wb, ...
      origconfigtrl, [], samprate_orig, struct(), samprate_orig, ...
      wave_plots_wanted, wave_zoom_ranges, wave_zoom_labels, inf, ...
      [ dataset.settitle ' - Waveforms - ' thistitle ], ...
      [ plotdir filesep dataset.setlabel '-waves-' thislabel ] );
  end

  if want_timelock_plots
    % Calculate the timelock average.
    ftdata_timelock = ft_timelockanalysis( struct(), ftdata_wb );

    euPlot_plotFTTimelock( ftdata_timelock, 2.0, timelock_plots_wanted, ...
      timelock_zoom_ranges_ms, timelock_zoom_labels, inf, ...
      [ dataset.settitle ' - Average Response - ' thistitle ], ...
      [ plotdir filesep dataset.setlabel '-avg-' thislabel ] );
  end


  %
  % Spectrogram plots.

  if want_spect_plots
    % Calculate per-trial spectrograms. This may take time!
    % FIXME - A constant-Q spectrogram would be more useful here.
    [ spectfreqs specttimes spectpowers ] = ...
      nlFT_calcSteppedWindowSpectrogram( ftdata_wb, ...
        spectrogram_win_secs, spectrogram_step_secs, ...
        [ min(ftdata_wb.time{1}) max(ftdata_wb.time{1}) ], ...
        spectrogram_freq_range );

    % Collapse this into an average spectrogram.
    % Summing is fine, since it's arbitrary units.

    freqcount = length(spectfreqs);
    timecount = length(specttimes);

    avgspect = sum(spectpowers, 1);
    avgspect = reshape(avgspect, chancount, freqcount, timecount);

    % FIXME - Move zoomed 2D plotting to a helper or library function at some point.

    thisfig = figure();
    figure(thisfig);

    for cidx = 1:chancount
      thisslice = avgspect(cidx,:,:);
      thisslice = reshape(thisslice, freqcount, timecount);

      for zidx = 1:length(spect_zoom_ranges)
        clf('reset');

        % FIXME - Auto-ranging axes.
        nlPlot_axesPlotSurface2D( gca, thisslice, specttimes, spectfreqs, ...
          spect_zoom_ranges{zidx}, spect_plot_freqrange, ...
          'linear', 'log', spect_loglin, 'Time (s)', 'Frequency (Hz)', ...
          [ dataset.settitle ' - Spectrogram - ' thistitle ...
            ' - ' chantitles{cidx} ' - ' spect_zoom_labels{zidx} ] );

        saveas( thisfig, ...
          [ plotdir filesep dataset.setlabel '-spect-' thislabel ...
            '-' chanlabels{cidx} '-' spect_zoom_labels{zidx} '.png' ] );
      end
    end

    close(thisfig);
  end


  %
  % Apply narrow-band filtering if desired.

  % This treats the signal as a LFP, even if we were reading MUA data.

  % Make a copy whether or not we're filtering it.
  ftdata_bpf = ftdata_wb;

  % FIXME - Get a nominal midband frequency even without the bandpass filter.
  % The Robinson/Hindriks stimulation has a peak at about 8 Hz.
  midband_freq = 8.0;

  if ~isempty(bandpass_band)

    midband_freq = geomean(bandpass_band);


    % Hindriks 2024 used a 4th-order Butterworth filter.
    % NOTE - FT complains about unstable filters. Specify a method to
    % stabilize them (the alternative would be to trap errors).
    filtconfig = struct( 'bpfilter', 'yes', 'bpfilttype', 'but', ...
      'bpfreq', [ min(bandpass_band), max(bandpass_band) ], ...
      'bpinstabilityfix', 'split', ...
      'bpfiltord', 4, 'feedback', 'no' );

    ftdata_bpf = ft_preprocessing( filtconfig, ftdata_wb );


    % Plot the filtered signals if we want waveform plots.

    if want_wave_plots
      euPlot_plotFTTrials( ftdata_bpf, samprate_wb, ...
        origconfigtrl, [], samprate_orig, struct(), samprate_orig, ...
        wave_plots_wanted, wave_zoom_ranges, wave_zoom_labels, inf, ...
        [ dataset.settitle ' - Narrow-Band - ' thistitle ], ...
        [ plotdir filesep dataset.setlabel '-lfp-' thislabel ] );
    end


    % FIXME - Cut and pasted spectrogram code. Make a helper!

    if want_spect_plots
      % Calculate per-trial spectrograms. This may take time!
      % FIXME - A constant-Q spectrogram would be more useful here.
      [ spectfreqs specttimes spectpowers ] = ...
        nlFT_calcSteppedWindowSpectrogram( ftdata_bpf, ...
          spectrogram_win_secs, spectrogram_step_secs, ...
          [ min(ftdata_bpf.time{1}) max(ftdata_bpf.time{1}) ], ...
          spectrogram_freq_range );

      % Collapse this into an average spectrogram.
      % Summing is fine, since it's arbitrary units.

      freqcount = length(spectfreqs);
      timecount = length(specttimes);

      avgspect = sum(spectpowers, 1);
      avgspect = reshape(avgspect, chancount, freqcount, timecount);

      % FIXME - Move zoomed 2D plotting to a helper or library function at some point.

      thisfig = figure();
      figure(thisfig);

      for cidx = 1:chancount
        thisslice = avgspect(cidx,:,:);
        thisslice = reshape(thisslice, freqcount, timecount);

        for zidx = 1:length(spect_zoom_ranges)
          clf('reset');

          % FIXME - Auto-ranging axes.
          nlPlot_axesPlotSurface2D( gca, thisslice, specttimes, spectfreqs, ...
            spect_zoom_ranges{zidx}, spect_plot_freqrange, ...
            'linear', 'log', spect_loglin, 'Time (s)', 'Frequency (Hz)', ...
            [ dataset.settitle ' - Spectrogram - ' thistitle ...
              ' - ' chantitles{cidx} ' - ' spect_zoom_labels{zidx} ] );

          saveas( thisfig, ...
            [ plotdir filesep dataset.setlabel '-lfpspect-' thislabel ...
              '-' chanlabels{cidx} '-' spect_zoom_labels{zidx} '.png' ] );
        end
      end

      close(thisfig);
    end
  end


  %
  % Phase statistics.

  % We want the statistics with or without a report, to tune parameters.

  phasedata = eiCalc_calcTrialPhaseStats( ftdata_bpf, ftdata_bpf, ...
    xcorr_params_phase, { 'pertrial' } );

  % Use a fixed threshold if provided or make an adaptive one if we have
  % a threshold. Otherwise accept all PLVs.

  if isnan(xcorr_phase_min_plv)
    xcorr_phase_min_plv = -1;

    if ~isnan(xcorr_phase_plv_prctile)
      [ plvmean plvdev plvprc ] = nlProc_calcMatrixStats( ...
        phasedata.plvtrials, [ xcorr_phase_plv_prctile ] );
      xcorr_phase_min_plv = plvprc;
    end
  end


  if want_phase_stats_report
    % Tattle the threshold here, too.
    disp(sprintf( '.. Using PLV threshold of %.2f.', xcorr_phase_min_plv ));

    % Get aggregate statistics.
    phasereporttext = doFT_reportPhaseStats( phasedata );

    % Add banners. Remember that the report itself already has \n.
    phasereporttext = sprintf( '%s\n%s%s\n', ...
      [ '-- Phase alignment report for case "' thislabel '":' ], ...
      phasereporttext,  '-- End of phase alignment report.' );

    nlIO_writeTextFile( [ plotdir filesep dataset.setlabel ...
      '-phase-stats-' thislabel '.txt' ], phasereporttext );

    disp(phasereporttext);

    globalphasestatreport = [ globalphasestatreport phasereporttext ];

  end



  %
  % Information flow plots.

  for infidx = 1:infocount

    thisinfolabel = infocases{infidx,1};
    thisinfotitle = infocases{infidx,2};
    thisinfoabbrv = infocases{infidx,3};
    thisinfobase = infocases{infidx,4};
    thisinfofunc = infocases{infidx,5};
    thisinfofuncphase = infocases{infidx,6};


    % FIXME - Kludge.
    thisinfoprefix = [ thisinfobase 'avg' ];
    if want_span_trials
      thisinfoprefix = [ thisinfobase 'concat' ];
    end

    thisinfodatafield = [ thisinfoprefix 'data' ];


    % Get info-specific plotting metadata.

    this_network_plotconfig = xcorr_network_config;

    if isfield( xcorr_network_amp_tuning, thisinfolabel )
      this_network_plotconfig.amp_range = ...
        xcorr_network_amp_tuning.(thisinfolabel);
    end


    if want_info_plots_lagtime || want_info_plots_line ...
      || want_info_plots_search || want_info_plots_stats_box ...
      || want_info_plots_stats_grid || want_info_plots_stats_report ...
      || want_info_plots_network || want_info_plots_network_anim ...
      || want_info_plots_byphase


      %
      % Calculate pairwise information and derived values.

      % Raw pariwise data data.

      disp([ '.. Computing ' thisinfotitle '.' ]);
      tic;

      xcorrdata = ...
        thisinfofunc( ftdata_bpf, ftdata_bpf, xcorr_params );

      durstring = nlUtil_makePrettyTime(toc);
      disp([ '.. Computed ' thisinfotitle ' in ' durstring '.' ]);

      if want_info_plots_byphase
        disp([ '.. Computing phase-binned ' thisinfotitle '.' ]);
        tic;

        [ xcorrphasetime xcorrphaselagboth xcorrphaselagafter ] = ...
        doFT_calcPhaseBinnedInfo( ...
          ftdata_bpf, ftdata_bpf, xcorr_params_phase, ...
          thisinfofuncphase, thisinfoprefix, xcorr_phase_min_plv, ...
          xcorr_phase_bins, xcorr_phase_bin_width_deg, ...
          xcorr_phase_time_smooth_ms, [ xcorr_after_start_ms inf ] );

        durstring = nlUtil_makePrettyTime(toc);
        disp([ '.. Computed phase-binned ' thisinfotitle ' in ' ...
          durstring '.' ]);

        % Diagnostics.
        disp(sprintf( '.. Used %d phase bins, %d degrees wide.', ...
          xcorr_phase_bins, xcorr_phase_bin_width_deg ));
      end


      % Processed pairwise information data without frequency or phase.

      % This returns struct arrays.
      % For vs lag, segment into "before" and "after" time series.
      [ xcorrvstime xcorrvslag ] = eiCalc_collapseTimeLagAverages( ...
        xcorrdata, thisinfoprefix, ...
        {[], [-inf xcorr_before_start_ms], [ xcorr_after_start_ms inf ]}, ...
        {[]} );

      % Vs time is a single element. Sort out vs lag.
      xcorrvslagbefore = xcorrvslag(2);
      xcorrvslagafter = xcorrvslag(3);
      xcorrvslag = xcorrvslag(1);

      % Get nominal peak locations.
      xcorr_peak_smooth_ms = ...
        xcorr_peak_smooth_wavecount * 1000 / midband_freq;

      xcorrpeaks = eiCalc_findTimeLagPeaks( xcorrdata, thisinfoprefix, ...
        xcorr_peak_smooth_ms, xcorr_peak_location, xcorr_peak_method );

      % Suppress peaks that weren't strong enough.
      if ~isnan(xcorr_time_amp_threshold)
        squashmask = ( abs(xcorrpeaks.peakamps) < xcorr_time_amp_threshold );
        xcorrpeaks.peakamps(squashmask) = NaN;
        xcorrpeaks.peaklags(squashmask) = NaN;
      end

      % Get nominal before/after peak statistics.

      [ xcorrampmeanbefore xcorrampdevbefore ...
        xcorrlagmeanbefore xcorrlagdevbefore ] = ...
        eiCalc_getTimeLagPeakStats( xcorrdata, thisinfoprefix, ...
          [ -inf xcorr_before_start_ms ], xcorr_peak_smooth_ms, ...
          xcorr_peak_halo_thresh, xcorr_peak_mag_accept );

      [ xcorrampmeanafter xcorrampdevafter ...
        xcorrlagmeanafter xcorrlagdevafter ] = ...
        eiCalc_getTimeLagPeakStats( xcorrdata, thisinfoprefix, ...
          [ xcorr_after_start_ms inf ], xcorr_peak_smooth_ms, ...
          xcorr_peak_halo_thresh, xcorr_peak_mag_accept );

% FIXME - Diagnostics.
if false
disp(sprintf( 'xx  Raw amp range:  %+.2e..%+.2e', ...
min(xcorrdata.(thisinfodatafield), [], 'all'), ...
max(xcorrdata.(thisinfodatafield), [], 'all') ));
disp(sprintf( 'xx  Mean amp before:  %+.2e..%+.2e   After:  %+.2e..%+.2e', ...
min(xcorrampmeanbefore, [], 'all'), max(xcorrampmeanbefore, [], 'all'), ...
min(xcorrampmeanafter, [], 'all'), max(xcorrampmeanafter, [], 'all') ));
end

      % Filter to get clean statistics.

      [ xcorrampgoodbefore xcorrlaggoodbefore ] = eiAux_makeValidMask( ...
        xcorrampmeanbefore, xcorrampdevbefore, ...
        xcorrlagmeanbefore, xcorrlagdevbefore, xcorr_stats_config );

      [ xcorrampgoodafter xcorrlaggoodafter ] = eiAux_makeValidMask( ...
        xcorrampmeanafter, xcorrampdevafter, ...
        xcorrlagmeanafter, xcorrlagdevafter, xcorr_stats_config );

      % Mask off "bad" cells.
      % We only need to do this for the mean, not deviation.

      if xcorr_squash_if_any_nan
        % Old behavior: If any statistic is bad for a cell, squash the cell.
        xcorrgoodmask = xcorrampgoodbefore & xcorrampgoodafter ...
          & xcorrlaggoodbefore & xcorrlaggoodafter;

        xcorrampgoodbefore = xcorrgoodmask;
        xcorrampgoodafter = xcorrgoodmask;
        xcorrlaggoodbefore = xcorrgoodmask;
        xcorrlaggoodafter = xcorrgoodmask;
      end

      xcorrampmeanbeforemasked = xcorrampmeanbefore;
      xcorrampmeanbeforemasked(~xcorrampgoodbefore) = NaN;

      xcorrampmeanaftermasked = xcorrampmeanafter;
      xcorrampmeanaftermasked(~xcorrampgoodafter) = NaN;

      xcorrlagmeanbeforemasked = xcorrlagmeanbefore;
      xcorrlagmeanbeforemasked(~xcorrlaggoodbefore) = NaN;

      xcorrlagmeanaftermasked = xcorrlagmeanafter;
      xcorrlagmeanaftermasked(~xcorrlaggoodafter) = NaN;

      % Get significant changes.

      xcorrampsigmask = nlProc_getSignificantChanges( ...
        xcorrampmeanbeforemasked, xcorrampdevbefore, ...
        xcorrampmeanaftermasked, xcorrampdevafter, ...
        xcorr_stats_config.significant_amp_sigma );

      xcorrlagsigmask = nlProc_getSignificantChanges( ...
        xcorrlagmeanbeforemasked, xcorrlagdevbefore, ...
        xcorrlagmeanaftermasked, xcorrlagdevafter, ...
        xcorr_stats_config.significant_lag_sigma );


      % If we were given an amplitude span of 'auto', auto-detect it.
      % FIXME - We should auto-range across all data cases, not just this one.
      if ischar(this_network_plotconfig.amp_range)
        % We have to squash the diagonal for XC or MI to make sense.
        % TE self-squashes that already.
        scratchbefore = xcorrampmeanbefore;
        scratchafter = xcorrampmeanafter;
        scratchcount = size(scratchbefore);
        scratchcount = min( scratchcount([ 1 2 ]) );
        for scratchidx = 1:scratchcount
          scratchbefore(scratchidx,scratchidx,:) = NaN;
          scratchafter(scratchidx,scratchidx,:) = NaN;
        end

        % Do the auto-ranging.
        scratch = reshape(scratchbefore, 1, []);
        scratch = [ scratch reshape(scratchafter, 1, []) ];
        scratch = scratch(~isnan(scratch));
        this_network_plotconfig.amp_range = ...
          prctile( scratch, xcorr_network_config.amp_auto_prctiles );

        disp(sprintf( ...
          '.. Auto-tuning network plot range for %s to %+.2e..%+.2e.', ...
          thisinfoabbrv, min(this_network_plotconfig.amp_range), ...
          max(this_network_plotconfig.amp_range) ));
      end


      %
      % Build various plotting labels and titles.
      % FIXME - Make a helper for "sprintf and prefix N and P for numbers".

      delaylist_ms = xcorrdata.delaylist_ms;
      delaycount = length(delaylist_ms);
      % FIXME - Blithely assume mostly-uniform spacing.
      delaystep_ms = median(diff( delaylist_ms ));

      [ delaylabels delaytitles ] = ...
        nlUtil_makeIntegerTimeLabels( delaylist_ms, 'ms' );

      timelist_ms = xcorrdata.windowlist_ms;
      timecount = length(timelist_ms);
      % FIXME - Blithely assume mostly-uniform spacing.
      timestep_ms = median(diff( timelist_ms ));

      [ timelabels timetitles ] = ...
        nlUtil_makeIntegerTimeLabels( timelist_ms, 'ms' );


      % Get the global Z data range, tolerating the "it's all zero" case.
      % NOTE - Ignore self-comparison.
      zrange = 0;
      for firstidx = 1:(chancount-1)
        for secondidx = (firstidx+1):chancount
          thisslice = ...
            xcorrdata.(thisinfodatafield)(firstidx,secondidx,:,:,:);
          thisslice = reshape( thisslice, timecount, delaycount );
          zrange = max( zrange, max( abs(thisslice), [], 'all' ) );
        end
      end
      zrange = max(zrange, 1e-3);
      zrange = [ -zrange zrange ];


      % There's no point in auto-ranging collapsed versions.
      % It would help if different test cases had different ranges, but
      % we're auto-ranging separately for every test case (set of trials).


      %
      % Render plots.


      % Common plotting.

      thisfig = figure();
      figure(thisfig);


      % Pairwise information vs lag and time for each channel pair.

      if want_info_plots_lagtime
        for firstidx = 1:(chancount-1)
          for secondidx = (firstidx+1):chancount
            % FIXME - Only plotting the first trial (we merged trials).
            trialidx = 1;

            thisslice = ...
              xcorrdata.(thisinfodatafield)(firstidx,secondidx,trialidx,:,:);
            thisslice = reshape( thisslice, timecount, delaycount );
            % Remember that we want (y,x).
            thisslice = transpose(thisslice);

            % Tolerate gaps, in case we're using this with real data.
            [ thisslice, thistimeseries, thislagseries ] = ...
              nlProc_padHeatmapGaps( thisslice, timelist_ms, delaylist_ms );

            clf('reset');
            colormap(thisfig, colmapamp);

            % FIXME - Auto-ranging axes.
            nlPlot_axesPlotSurface2D( gca, thisslice, ...
              thistimeseries, thislagseries, [], [], ...
              'linear', 'linear', 'linear', ...
              'Time (ms)', 'Delay (ms)', ...
              [ dataset.settitle ' - ' thisinfoabbrv ' - ' thistitle ...
              ' - ' chantitles{firstidx} ' and ' chantitles{secondidx} ] );

            clim(zrange);

            thiscol = colorbar;
            thiscol.Label.String = [ 'Normalized ' thisinfotitle ];

            saveas( thisfig, [ plotdir filesep dataset.setlabel ...
                '-' thisinfolabel '-lagtime-' thislabel ...
                '-' chanlabels{firstidx} '-' chanlabels{secondidx} '.png' ] );
          end
        end
      end


      % Pairwise info vs lag _or_ time, collapsing the other, for each pair.

      if want_info_plots_line

        % Convert ms to samples. May be NaN, for no smoothing.
        smoothwindow = round( xcorr_collapse_smooth_ms / timestep_ms );


        % Vs window time.
        % These are indexed by (first,second,trialidx,winidx).

        doFT_plotPairWaves( xcorrvstime.dev, chantitles, chantitles, ...
          xcorrvstime.trialnums, ...
          timelist_ms, smoothwindow, [], xcorr_plot_devcursors, ...
          'Time (ms)', [ thisinfotitle ' Standard Deviation' ], ...
          [ dataset.settitle ' - ' thisinfoabbrv ...
            ' Dev vs Time - ' thistitle ], ...
          [ plotdir filesep dataset.setlabel ...
            '-' thisinfolabel '-vstimedev-' thislabel ] );

        if xcorr_plot_want_timemean
          doFT_plotPairWaves( xcorrvstime.avg, chantitles, chantitles, ...
            xcorrvstime.trialnums, ...
            timelist_ms, smoothwindow, [], xcorr_plot_meancursors, ...
            'Time (ms)', [ 'Mean ' thisinfotitle ], ...
            [ dataset.settitle ' - ' thisinfoabbrv ...
              ' Mean vs Time - ' thistitle ], ...
            [ plotdir filesep dataset.setlabel ...
              '-' thisinfolabel '-vstimeavg-' thislabel ] );
        end


        % Vs time lag.
        % These are indexed by (first,second,trialidx,lagidx).

        % NOTE - We have "before", "after", and "all" as time options.
        xcorrvslaglut = struct();
        xcorrvslaglut.before = xcorrvslagbefore;
        xcorrvslaglut.after = xcorrvslagafter;
        xcorrvslaglut.all = xcorrvslag;

        for pidx = 1:length(xcorr_plot_list_vslag)
          thiskey = xcorr_plot_list_vslag{pidx};
          thisdata = xcorrvslaglut.(thiskey);

          doFT_plotPairWaves( thisdata.avg, chantitles, chantitles, ...
            thisdata.trialnums, ...
            delaylist_ms, NaN, [], xcorr_plot_meancursors, ...
            'Delay (ms)', [ 'Mean ' thisinfotitle ], ...
            [ dataset.settitle ' - ' thisinfoabbrv ...
              ' Mean vs Lag (' thiskey ') - ' thistitle ], ...
            [ plotdir filesep dataset.setlabel ...
              '-' thisinfolabel '-vslagavg-' thiskey '-' thislabel ] );

          if xcorr_plot_want_lagdev
            doFT_plotPairWaves( thisdata.dev, chantitles, chantitles, ...
              thisdata.trialnums, ...
              delaylist_ms, NaN, [], xcorr_plot_devcursors, ...
              'Delay (ms)', [ thisinfotitle ' Standard Deviation' ], ...
              [ dataset.settitle ' - ' thisinfoabbrv ...
                ' Dev vs Lag (' thiskey ') - ' thistitle ], ...
              [ plotdir filesep dataset.setlabel '-' ...
                thisinfolabel '-vslagdev-' thiskey '-' thislabel ] );
          end
        end

      end


      % Pairwise information peak evolution over time, for each pair.

      if want_info_plots_search

        doFT_plotPairWaves( xcorrpeaks.peakamps, chantitles, chantitles, ...
          xcorrpeaks.trialnums, ...
          timelist_ms, NaN, [], xcorr_plot_ampcursors, ...
          'Time (ms)', [ 'Strongest ' thisinfotitle ], ...
          [ dataset.settitle ' - ' thisinfoabbrv ...
            ' Peak Amp vs Time - ' thistitle ], ...
          [ plotdir filesep dataset.setlabel ...
            '-' thisinfolabel '-peakamp-' thislabel ] );

        doFT_plotPairWaves( xcorrpeaks.peaklags, chantitles, chantitles, ...
          xcorrpeaks.trialnums, ...
          timelist_ms, NaN, [], xcorr_plot_lagcursors, ...
          'Time (ms)', 'Delay at Peak (ms)', ...
          [ dataset.settitle ' - ' thisinfoabbrv ...
            ' Peak Lag vs Time - ' thistitle ], ...
          [ plotdir filesep dataset.setlabel ...
            '-' thisinfolabel '-peaklag-' thislabel ] );

      end


      % Before-and-after statistics for each pair - box chart.

      if want_info_plots_stats_box
        % FIXME - Only plotting the first trial (we merged trials).
        trialidx = 1;

        confsigma = xcorr_stats_config.significant_amp_sigma;
        doFT_plotBeforeAfterStats( ...
          xcorrampmeanbefore(:,:,trialidx), ...
          confsigma * xcorrampdevbefore(:,:,trialidx), ...
          xcorrampmeanafter(:,:,trialidx), ...
          confsigma * xcorrampdevafter(:,:,trialidx), ...
          chantitles, chantitles, ...
          [], xcorr_plot_ampcursors, [ thisinfotitle ' Peak Amplitude' ], ...
          [ dataset.settitle ' - Average ' thisinfoabbrv ...
            ' Peak Amp - ' thistitle ], ...
          [ plotdir filesep dataset.setlabel ...
            '-' thisinfolabel '-statsamp-' thislabel '.png' ] );

        confsigma = xcorr_stats_config.significant_lag_sigma;
        doFT_plotBeforeAfterStats( ...
          xcorrlagmeanbefore(:,:,trialidx), ...
          confsigma * xcorrlagdevbefore(:,:,trialidx), ...
          xcorrlagmeanafter(:,:,trialidx), ...
          confsigma * xcorrlagdevafter(:,:,trialidx), ...
          chantitles, chantitles, ...
          [], xcorr_plot_lagcursors, [ thisinfotitle ' Peak Delay (ms)' ], ...
          [ dataset.settitle ' - Average ' thisinfoabbrv ...
            ' Peak Lag - ' thistitle ], ...
          [ plotdir filesep dataset.setlabel ...
            '-' thisinfolabel '-statslag-' thislabel '.png' ] );
      end


      % Filtered before-and-after statistics for each pair - report.

      if want_info_plots_stats_report

        xcorrreporttext = doFT_reportInfoStats( ...
          xcorrampmeanbeforemasked, xcorrampmeanaftermasked, ...
          xcorrlagmeanbeforemasked, xcorrlagmeanaftermasked, ...
          xcorrampsigmask, xcorrlagsigmask, chantitles, chantitles );

        % Add banners. Remember that the report itself already has \n.
        xcorrreporttext = sprintf( '%s\n%s%s\n', ...
          [ '-- ' thisinfotitle ' report for case "' thislabel '":' ], ...
          xcorrreporttext, [ '-- End of ' thisinfotitle ' report.' ] );

        nlIO_writeTextFile( [ plotdir filesep dataset.setlabel ...
          '-' thisinfolabel '-stats-' thislabel '.txt' ], ...
          xcorrreporttext );

        disp(xcorrreporttext);

        globalxcorrstatreport = [ globalxcorrstatreport xcorrreporttext ];

      end


      % Filtered before-and-after statistics for each pair - grid.

      if want_info_plots_stats_grid
        % FIXME - Only plotting the first trial (we merged trials).
        trialidx = 1;

        doFT_plotBeforeAfterGrid( ...
          xcorrampmeanbeforemasked(:,:,trialidx), ...
          xcorrampmeanaftermasked(:,:,trialidx), ...
          xcorrlagmeanbeforemasked(:,:,trialidx), ...
          xcorrlagmeanaftermasked(:,:,trialidx), ...
          xcorrampsigmask(:,:,trialidx), ...
          xcorrlagsigmask(:,:,trialidx), ...
          chantitles, chantitles, colmapamp, colmaplag, ...
          thisinfotitle, thisinfoabbrv, ...
          [ dataset.settitle ' - Average ' thisinfoabbrv ' - ' thistitle ], ...
          [ plotdir filesep dataset.setlabel ...
            '-' thisinfolabel '-statgrid-' thislabel '.png' ] );
      end


      % Network graph of the average pairwise information statistics.

      if want_info_plots_network
        % FIXME - Only plotting the first trial (we merged trials).
        trialidx = 1;

        euPlot_hlevPlotNetwork( ...
          xcorrampmeanbefore(:,:,trialidx), ...
          xcorrlagmeanbefore(:,:,trialidx), ...
          [], this_network_plotconfig, chantitles, chantitles, ...
          [ dataset.settitle ' - ' thisinfoabbrv ...
            ' Network Before - ' thistitle ], ...
          [ plotdir filesep dataset.setlabel ...
            '-' thisinfolabel '-netbefore-' thislabel '.png' ] );

        euPlot_hlevPlotNetwork( ...
          xcorrampmeanafter(:,:,trialidx), ...
          xcorrlagmeanafter(:,:,trialidx), ...
          [], this_network_plotconfig, chantitles, chantitles, ...
          [ dataset.settitle ' - ' thisinfoabbrv ...
            ' Network After - ' thistitle ], ...
          [ plotdir filesep dataset.setlabel ...
            '-' thisinfolabel '-netafter-' thislabel '.png' ] );
      end


      % Network graph showing how the detected peaks evolve over time.

      if want_info_plots_network_anim

        % FIXME - Only plotting the first trial (we merged trials).
        trialidx = 1;

        % FIXME - Hope that element order stays consistent!

        scratchamps = xcorrpeaks.peakamps(:,:,trialidx,:);
        scratchlags = xcorrpeaks.peaklags(:,:,trialidx,:);

        scratchsize = size(scratchamps);
        scratchsize = scratchsize([ 1 2 4 ]);

        scratchamps = reshape(scratchamps, scratchsize);
        scratchlags = reshape(scratchlags, scratchsize);


        euPlot_hlevPlotNetwork( ...
          scratchamps, scratchlags, timelist_ms, ...
          this_network_plotconfig, chantitles, chantitles, ...
          [ dataset.settitle ' - ' thisinfoabbrv ' Network - ' thistitle ], ...
          [ plotdir filesep dataset.setlabel ...
            '-' thisinfolabel '-netvideo-' thislabel '.avi' ] );

      end


      % Plots of pairwise info vs lag and phase difference, for each pair.

      if want_info_plots_byphase
        % NOTE - We have "full time range" and "after-stim only" ranges.

        if want_info_plots_byphase_more
          doFT_plotPairVsPhase( xcorrphaselagboth, xcorr_plot_phasenan, ...
            'avg', [ 'Average ' thisinfotitle ], ...
            'delaylist_ms', 'Delay (ms)', thisinfotitle, thisinfoabbrv, ...
            colmapamp, { 'heat', 'line' }, ...
            [ dataset.settitle ' - Avg ' thisinfoabbrv ...
              ' vs Phase Diff - ' thistitle ], ...
            [ plotdir filesep dataset.setlabel ...
              '-' thisinfolabel '-lagphaseboth-' thislabel ] );
        end

        doFT_plotPairVsPhase( xcorrphaselagafter, xcorr_plot_phasenan, ...
          'avg', [ 'Average ' thisinfotitle ], ...
          'delaylist_ms', 'Delay (ms)', thisinfotitle, thisinfoabbrv, ...
          colmapamp, { 'heat', 'line' }, ...
          [ dataset.settitle ' - Avg ' thisinfoabbrv ...
            ' vs Phase Diff after - ' thistitle ], ...
          [ plotdir filesep dataset.setlabel ...
            '-' thisinfolabel '-lagphaseafter-' thislabel ] );

        % NOTE - For averaging across lags, use RMS magnitude (deviation),
        % not the average (which should be zero).
        % In practice, neither is particularly informative.

        if want_info_plots_byphase_more
          doFT_plotPairVsPhase( xcorrphasetime, xcorr_plot_phasenan, ...
            'avg', [ 'Average ' thisinfotitle ], ...
            'windowlist_ms', 'Time (ms)', thisinfotitle, thisinfoabbrv, ...
            colmapamp, { 'heat' }, ...
            [ dataset.settitle ' - Avg ' thisinfoabbrv ...
              ' vs Phase Diff - ' thistitle ], ...
            [ plotdir filesep dataset.setlabel ...
              '-' thisinfolabel '-timephasemean-' thislabel ] );
        end

        % NOTE - For averaging across lags, use RMS magnitude (deviation),
        % not the average (which should be zero).
        doFT_plotPairVsPhase( xcorrphasetime, xcorr_plot_phasenan, ...
          'dev', [ 'RMS Average ' thisinfotitle ], ...
          'windowlist_ms', 'Time (ms)', thisinfotitle, thisinfoabbrv, ...
          colmapamp, { 'heat' }, ...
          [ dataset.settitle ' - RMS ' thisinfoabbrv ...
            ' vs Phase Diff - ' thistitle ], ...
          [ plotdir filesep dataset.setlabel ...
            '-' thisinfolabel '-timephasedev-' thislabel ] );
      end


      % Common plotting.

      close(thisfig);

    end


    %
    % Pairwise information vs frequency plots.

    if want_info_plots_byfreq

      % First pass: Aggregate data.

      timelist_ms = xcorr_params.timelist_ms;
      timecount = length(timelist_ms);
      bandcount = length(xcorr_vsfreq_bands);

      % Diagnostics.
      disp(sprintf( '.. Using %d frequency bins.', bandcount ));


      xcvsfreq_rms = nan([ chancount chancount timecount bandcount ]);
      xcvsfreq_peak = nan(size(xcvsfreq_rms));
      freqlist = nan(size(xcorr_vsfreq_bands));

      for bidx = 1:bandcount

        thisband = xcorr_vsfreq_bands{bidx};

        % Use the geometric mean, so that log scale bin edges reconstruct
        % properly.
        thismidfreq = geomean(thisband);
        freqlist(bidx) = thismidfreq;


        disp(sprintf( ...
          [ '.. Computing ' thisinfotitle ' for %d - %d Hz.' ], ...
          round(min(thisband)), round(max(thisband)) ));
        tic;


        % Filter the data.

        % Hindriks 2024 used a 4th-order Butterworth filter.
        % NOTE - FT complains about unstable filters. Specify a method to
        % stabilize them (the alternative would be to trap errors).
        filtconfig = struct( 'bpfilter', 'yes', 'bpfilttype', 'but', ...
          'bpfreq', [ min(thisband), max(thisband) ], ...
          'bpinstabilityfix', 'split', ...
          'bpfiltord', 4, 'feedback', 'no' );

        ftdata_nb = ft_preprocessing( filtconfig, ftdata_wb );


        % Get raw pairwise information.

        bandxcorrdata = ...
          thisinfofunc( ftdata_nb, ftdata_nb, xcorr_params );


        % Collapse this into a time series.
        % NOTE - Use _both_ RMS average _and_ peak detection, separately.

        % This returns struct arrays.
        % Use the entire time span when figuring out where to look for peaks
        % in the lag domain.
        [ bandxcorrvstime bandxcorrvslag ] = ...
          eiCalc_collapseTimeLagAverages( ...
            bandxcorrdata, thisinfoprefix, {[]}, {[]} );

        % Get nominal peak locations.
        xcorr_peak_smooth_ms = ...
          xcorr_peak_smooth_wavecount * 1000 / thismidfreq;
        bandxcorrpeaks = ...
          eiCalc_findTimeLagPeaks( bandxcorrdata, thisinfoprefix, ...
            xcorr_peak_smooth_ms, xcorr_peak_location, xcorr_peak_method );

        % Suppress peaks that weren't strong enough.
        if ~isnan(xcorr_time_amp_threshold)
          squashmask = ...
            ( abs(bandxcorrpeaks.peakamps) < xcorr_time_amp_threshold );
          bandxcorrpeaks.peakamps(squashmask) = NaN;
          bandxcorrpeaks.peaklags(squashmask) = NaN;
        end


        % Store this slice's data.
        % NOTE - This can get big.


        % Indexing is (first,second,trial,window). We want to drop trialidx.
        thisrms = bandxcorrvstime.dev;
        thispeak = bandxcorrpeaks.peakamps;


        % FIXME - Only plotting the first trial (we merged trials).
        % FIXME - Hope that element order stays consistent!

        trialidx = 1;

        thisrms = thisrms(:,:,trialidx,:);
        thispeak = thispeak(:,:,trialidx,:);

        scratchsize = size(thisrms);
        scratchsize = scratchsize([ 1 2 4 ]);

        thisrms = reshape(thisrms, scratchsize);
        thispeak = reshape(thispeak, scratchsize);


        % Indexing is (first,second,window,freq).
        xcvsfreq_rms(1:chancount,1:chancount,1:timecount,bidx) = thisrms;
        xcvsfreq_peak(1:chancount,1:chancount,1:timecount,bidx) = thispeak;

        durstring = nlUtil_makePrettyTime(toc);
        disp([ '.. Computed narrow-band ' thisinfotitle ' in ' ...
          durstring '.' ]);

      end


      %
      % Second pass: Plot XC vs frequency and time.

      zrange_rms = 0;
      zrange_peak = 0;

      pairmask = nlUtil_getPairMask(chanlabels, chanlabels);

      for bidx = 1:bandcount
        for widx = 1:timecount
          % Indexing is (first,second,window,freq).

          thisslice = xcvsfreq_rms(:,:,widx,bidx);
          thisslice(~pairmask) = NaN;
          scratch = max( abs(thisslice), [], 'all' );
          zrange_rms = max(zrange_rms, scratch);

          thisslice = xcvsfreq_peak(:,:,widx,bidx);
          thisslice(~pairmask) = NaN;
          scratch = max( abs(thisslice), [], 'all' );
          zrange_peak = max(zrange_peak, scratch);
        end
      end

      thisfig = figure();
      figure(thisfig);

      for firstidx = 1:(chancount-1)
        for secondidx = (firstidx+1):chancount

          % Indexing is (first,second,window,freq).
          % Remember that we want (y,x) out.

          thisrms = xcvsfreq_rms(firstidx,secondidx,:,:);
          thisrms = reshape( thisrms, timecount, bandcount);
          thisrms = transpose(thisrms);

          thispeak = xcvsfreq_peak(firstidx,secondidx,:,:);
          thispeak = reshape( thispeak, timecount, bandcount);
          thispeak = transpose(thispeak);


          clf('reset');
          colormap(thisfig, colmapamp);

          % Tolerate gaps in the time axis, in case we're using real data.
          % NOTE - We have to make the frequency axis linear to call this.
          [ thisdata, thistimeseries, thisfreqseries ] = ...
            nlProc_padHeatmapGaps( thisrms, timelist_ms, log(freqlist));
          thisfreqseries = exp(thisfreqseries);

          nlPlot_axesPlotSurface2D( gca, thisdata, ...
            thistimeseries, thisfreqseries, [], [], ...
            'linear', 'log', 'linear', ...
            'Time (ms)', 'Frequency (Hz)', ...
            [ dataset.settitle ' - RMS ' thisinfoabbrv ' vs Freq - ' ...
              thistitle ' - ' chantitles{firstidx} ' and ' ...
              chantitles{secondidx} ] );

          clim([ -zrange_rms, zrange_rms ]);

          thiscol = colorbar;
          thiscol.Label.String = [ 'RMS ' thisinfotitle ];

          saveas( thisfig, [ plotdir filesep dataset.setlabel ...
            '-' thisinfolabel '-freqtime-rms-' thislabel ...
            '-' chanlabels{firstidx} '-' chanlabels{secondidx} '.png' ] );


          clf('reset');
          colormap(thisfig, colmapamp);

          % Tolerate gaps in the time axis, in case we're using real data.
          % NOTE - We have to make the frequency axis linear to call this.
          [ thisdata, thistimeseries, thisfreqseries ] = ...
            nlProc_padHeatmapGaps( thispeak, timelist_ms, log(freqlist));
          thisfreqseries = exp(thisfreqseries);

          nlPlot_axesPlotSurface2D( gca, thisdata, ...
            thistimeseries, thisfreqseries, [], [], ...
            'linear', 'log', 'linear', ...
            'Time (ms)', 'Frequency (Hz)', ...
            [ dataset.settitle ' - Peak ' thisinfoabbrv ' vs Freq - ' ...
              thistitle ' - ' chantitles{firstidx} ' and ' ...
              chantitles{secondidx} ] );

          clim([ -zrange_peak, zrange_peak ]);

          thiscol = colorbar;
          thiscol.Label.String = [ 'Peak ' thisinfotitle ];

          saveas( thisfig, [ plotdir filesep dataset.setlabel ...
            '-' thisinfolabel '-freqtime-peak-' thislabel ...
            '-' chanlabels{firstidx} '-' chanlabels{secondidx} '.png' ] );

        end
      end

      close(thisfig);

    end

  end


  %
  % Hindriks 2023 analyses - Coherence, power correlation, and co-kurtosis.

  if want_htplots_coherence || want_htplots_powercorrel ...
    || want_htplots_powernongauss

    % This squashes the time delay range (only one lag value).
    % So, data matrices are (firstchans,secondchans,1,wintimes,1).

    htdata = eiCalc_calcHTStats( ...
      ftdata_bpf, ftdata_bpf, xcorr_params, { 'avgtrials' } );


    timelist_ms = htdata.windowlist_ms;
    timecount = length(timelist_ms);

    [ timelabels timetitles ] = ...
      nlUtil_makeIntegerTimeLabels( timelist_ms, 'ms' );


    % Get the global data ranges, tolerating the "data is all zero" case.
    % NOTE - Ignore self-comparison.

    zrange_coherence = 0;
    zrange_powercorrel = 0;
    zrange_nongauss = 0;

    % FIXME - Only plotting the first trial (we merged trials).
    trialidx = 1;

    for firstidx = 1:(chancount-1)
      for secondidx = (firstidx+1):chancount

        thisslice = htdata.coherenceavgdata(firstidx,secondidx,trialidx,:,:);
        thisslice = reshape( thisslice, 1, timecount );
        zrange_coherence = max( zrange_coherence, max(abs(thisslice)) );

        thisslice = htdata.powercorrelavgdata(firstidx,secondidx,trialidx,:,:);
        thisslice = reshape( thisslice, 1, timecount );
        zrange_powercorrel = max( zrange_powercorrel, max(abs(thisslice)) );

        thisslice = htdata.nongaussavgdata(firstidx,secondidx,trialidx,:,:);
        thisslice = reshape( thisslice, 1, timecount );
        zrange_nongauss = max( zrange_nongauss, max(abs(thisslice)) );

      end
    end

    zrange_coherence = max(zrange_coherence, 1e-3);
    zrange_powercorrel = max(zrange_powercorrel, 1e-3);
    zrange_nongauss = max(zrange_nongauss, 1e-3);

    zrange_coherence = [ -zrange_coherence zrange_coherence ];
    zrange_powercorrel = [ -zrange_powercorrel zrange_powercorrel ];
    zrange_nongauss = [ -zrange_nongauss zrange_nongauss ];


    % Statistics vs time for each channel pair.
    % One plot per first channel; second channel x time in each plot.

    % FIXME - Move 2D plotting to a helper or library function at some point.

    thisfig = figure();
    figure(thisfig);

    chanindices = 1:chancount;

    % FIXME - Only plotting the first trial (we merged trials).
    trialidx = 1;


    if want_htplots_coherence
      for firstidx = 1:chancount
        thisslice = htdata.coherenceavgdata(firstidx,:,trialidx,:,:);
        thisslice = reshape( thisslice, chancount, timecount );

        % We're already in (y,x) form; no need to transpose.

        % Squash self-comparison.
        thisslice(firstidx,:) = NaN;

        % NOTE - Coherence is complex; turn it into a magnitude.
        thisslice = abs(thisslice);

        % Tolerate gaps, in case we're using this with real data.
        [ thisslice, thistimeseries, thischanseries ] = ...
          nlProc_padHeatmapGaps( thisslice, timelist_ms, chanindices );

        clf('reset');
        colormap(thisfig, colmapamp);

        % FIXME - Auto-ranging axes.
        nlPlot_axesPlotSurface2D( gca, thisslice, ...
          thistimeseries, chanlabels, [], [], ...
          'linear', 'linear', 'linear', ...
          'Time (ms)', 'Channel', ...
          [ dataset.settitle ' - HT Coherence - ' thistitle ...
          ' - ' chantitles{firstidx} ] );

        % Flip the vertical axis so labels are in lexical order.
        set(gca, 'YDir', 'reverse');

        clim(zrange_coherence);

        thiscol = colorbar;
        thiscol.Label.String = 'abs( Coherence )';

        saveas( thisfig, ...
          [ plotdir filesep dataset.setlabel '-htcoherence-' thislabel ...
            '-' chanlabels{firstidx} '.png' ] );
      end
    end


    if want_htplots_powercorrel
      for firstidx = 1:chancount
        thisslice = htdata.powercorrelavgdata(firstidx,:,trialidx,:,:);
        thisslice = reshape( thisslice, chancount, timecount );

        % We're already in (y,x) form; no need to transpose.

        % Squash self-comparison.
        thisslice(firstidx,:) = NaN;

        % Tolerate gaps, in case we're using this with real data.
        [ thisslice, thistimeseries, thischanseries ] = ...
          nlProc_padHeatmapGaps( thisslice, timelist_ms, chanindices );

        clf('reset');
        colormap(thisfig, colmapamp);

        % FIXME - Auto-ranging axes.
        nlPlot_axesPlotSurface2D( gca, thisslice, ...
          thistimeseries, chanlabels, [], [], ...
          'linear', 'linear', 'linear', ...
          'Time (ms)', 'Channel', ...
          [ dataset.settitle ' - HT Power Correlation - ' thistitle ...
          ' - ' chantitles{firstidx} ] );

        % Flip the vertical axis so labels are in lexical order.
        set(gca, 'YDir', 'reverse');

        clim(zrange_powercorrel);

        thiscol = colorbar;
        thiscol.Label.String = 'Power Correlation';

        saveas( thisfig, ...
          [ plotdir filesep dataset.setlabel '-htpowercorrel-' thislabel ...
            '-' chanlabels{firstidx} '.png' ] );
      end
    end


    if want_htplots_powernongauss
      for firstidx = 1:chancount
        thisslice = htdata.nongaussavgdata(firstidx,:,trialidx,:,:);
        thisslice = reshape( thisslice, chancount, timecount );

        % We're already in (y,x) form; no need to transpose.

        % Squash self-comparison.
        thisslice(firstidx,:) = NaN;

        % Tolerate gaps, in case we're using this with real data.
        [ thisslice, thistimeseries, thischanseries ] = ...
          nlProc_padHeatmapGaps( thisslice, timelist_ms, chanindices );

        clf('reset');
        colormap(thisfig, colmapamp);

        % FIXME - Auto-ranging axes.
        nlPlot_axesPlotSurface2D( gca, thisslice, ...
          thistimeseries, chanlabels, [], [], ...
          'linear', 'linear', 'linear', ...
          'Time (ms)', 'Channel', ...
          [ dataset.settitle ' - HT Non-Gaussian Correl - ' thistitle ...
          ' - ' chantitles{firstidx} ] );

        % Flip the vertical axis so labels are in lexical order.
        set(gca, 'YDir', 'reverse');

        clim(zrange_nongauss);

        thiscol = colorbar;
        thiscol.Label.String = 'NG Power Correl';

        saveas( thisfig, ...
          [ plotdir filesep dataset.setlabel '-htnongauss-' thislabel ...
            '-' chanlabels{firstidx} '.png' ] );
      end
    end


    close(thisfig);

  end


  % Done.

end

disp([ '== Finished processing "' dataset.settitle '".' ]);



% Save and/or plot any global data.

if ~isempty(globalxcorrstatreport)
  nlIO_writeTextFile( ...
    [ plotdir filesep dataset.setlabel '-xc-stats.txt' ], ...
    globalxcorrstatreport );
end

if ~isempty(globalphasestatreport)
  nlIO_writeTextFile( ...
    [ plotdir filesep dataset.setlabel '-phase-stats.txt' ], ...
    globalphasestatreport );
end



%
% This is the end of the file.
