function doFT_plotPairVsPhase( xcorrdata, nanfill, datafield, datatitle, ...
  timefield, timetitle, infotitle, infoabbrv, ...
  colmap, plotswanted, titleprefix, fnameprefix )

% function doFT_plotPairVsPhase( xcorrdata, nanfill, datafield, datatitle, ...
%   timefield, timetitle, infotitle, infoabbrv, ...
%   colmap, plotswanted, titleprefix, fnameprefix )
%
% This plots pairwise information amplitude or deviation between each pair of
% channels as a function of time or lag and the phase difference between the
% channels' signals.
%
% This optionally also collapses the time/lag axis to make a line plot vs
% phase.
%
% "xcorrdata" is a structure produced by doFT_calcPhaseBinnedInfo. Data
%   matrices are indexed by (destchan, srcchan, trial, timeindex, phaseindex).
% "nanfill" is a value to replace NaN with in the data structure. This may
%   itself be NaN.
% "datafield" is the field to look for data values in (typically 'xcorravg'
%   or 'xcorrdev').
% "datatitle" is the title to use for the colour (Z) axis.
% "timefield" is the field to look for X axis coordinates in (typically
%   'delaylist_ms' or 'windowlist_ms').
% "timetitle" is the title to use for the X axis.
% "infotitle" is the name of the data value being plotted (e.g.
%   'Cross-Correlation').
% "infoabbrv" is an abbreviation of the name of the data value (e.g. 'XC').
% "colmap" is the colour map to use when plotting.
% "plotswanted" is a cell array containing zero or more of 'line' and 'heat'.
% "titleprefix" is a character vector used for building plot titles.
% "fnameprefix" is a character vector used for building plot filenames.
%
% No return value.


%
% Get unpacked data.

destlabels = xcorrdata.destchans;
srclabels = xcorrdata.srcchans;

[ destlabels desttitles ] = euUtil_makeSafeStringArray(destlabels);
[ srclabels srctitles ] = euUtil_makeSafeStringArray(srclabels);

phasevals = xcorrdata.phaselist_deg;
timevals = xcorrdata.(timefield);

trialnums = xcorrdata.trialnums;

datavals = xcorrdata.(datafield);

% Swap out NaN.
datavals(isnan(datavals)) = nanfill;


%
% Get metadata.

destcount = length(destlabels);
srccount = length(srclabels);

phasecount = length(phasevals);
timecount = length(timevals);

trialcount = size(xcorrdata.avg,3);

pairmask = nlUtil_getPairMask(destlabels, srclabels);

want_heat = any(strcmp( plotswanted, 'heat' ));
want_line = any(strcmp( plotswanted, 'line' ));


%
% Get the data range.

% Make sure to mask off the diagonal first.
for destidx = 1:destcount
  for srcidx = 1:srccount
    if ~pairmask(destidx,srcidx)
      datavals(destidx,srcidx,:,:) = NaN;
    end
  end
end

% This tolerates NaN data.
scratch = max( abs(datavals), [], 'all' );
scratch = max(scratch, 1e-30);
datarange = [ (-scratch) scratch ];


%
% Generate plots.

thisfig = figure();
figure(thisfig);

cols = nlPlot_getColorPalette();

for destidx = 1:destcount
  for srcidx = 1:srccount
    if pairmask(destidx,srcidx)
      for trialidx = 1:trialcount

        % Get heatmap data.

        thisslice = datavals(destidx,srcidx,trialidx,:,:);
        thisslice = reshape( thisslice, timecount, phasecount );
        % We want (y,x). It's already set up that way; X is phase, Y time.

        % Collapse this to get line plot data.

        rmsvals = [];
        peakvals = [];
        for pidx = 1:phasecount
          thiscol = thisslice(:,pidx);
          rmsvals(pidx) = std(thiscol);

          % FIXME - Kludge peak detection. No tracking or smoothing or window.
          thismin = min(thiscol);
          thismax = max(thiscol);
          if abs(thismin) > abs(thismax)
            peakvals(pidx) = thismin;
          else
            peakvals(pidx) = thismax;
          end
        end


        % Build trial suffixes.

        trialtitle = '';
        triallabel = '';

        if trialcount > 1
          trialtitle = sprintf( ' - Tr %04d', trialnums(trialidx) );
          triallabel = sprintf( '-tr%04d', trialnums(trialidx) );
        end


        % Heatmap plot.

        if want_heat
          % Tolerate discontinuous axes.
          [ thisslice, thisphaseseries, thistimeseries ] = ...
            nlProc_padHeatmapGaps( thisslice, phasevals, timevals );

          % Don't worry about NaNs. The plotting function handles them.

          clf('reset');
          colormap(thisfig, colmap);

          nlPlot_axesPlotSurface2D( gca, thisslice, ...
            thisphaseseries, thistimeseries, [], [], ...
            'linear', 'linear', 'linear', 'Phase (degrees)', timetitle, ...
            [ titleprefix ' - ' desttitles{destidx} ' from ' ...
              srctitles{srcidx} trialtitle ] );

          clim(datarange);

          thiscol = colorbar;
          thiscol.Label.String = datatitle;

          saveas( thisfig, [ fnameprefix '-' destlabels{destidx} ...
            '-' srclabels{srcidx} triallabel '.png' ] );
        end


      % Line plot.

      if want_line
        clf('reset');

        hold on;

        plot( phasevals, peakvals, 'Color', cols.blu, ...
          'DisplayName', [ 'Peak ' infoabbrv ] );
        plot( phasevals, rmsvals, 'Color', cols.brn, ...
          'DisplayName', [ 'RMS Avg ' infoabbrv ] );

        plot( phasevals, zeros(size(phasevals)), 'Color', cols.gry, ...
          'HandleVisibility', 'off' );

        hold off;

        ylim(datarange);

        xlabel('Phase (degrees)');
        ylabel(infotitle);

        legend('Location', 'southeast');

        title([ titleprefix ' - ' desttitles{destidx} ' from ' ...
          srctitles{srcidx} trialtitle ]);

        saveas( thisfig, [ fnameprefix '-line-' destlabels{destidx} ...
          '-' srclabels{srcidx} triallabel '.png' ] );
      end

      end
    end
  end
end


% Done.
end


%
% This is the end of the file.
