function doFT_plotBeforeAfterStats( ...
  meanbefore, confbefore, meanafter, confafter, ...
  firstlabels, secondlabels, yrange, ycursors, ytitle, figtitle, fname )

% function doFT_plotBeforeAfterStats( ...
%   meanbefore, confbefore, meanafter, confafter, ...
%   firstlabels, secondlabels, yrange, ycursors, ytitle, figtitle, fname )
%
% This makes a plot showing how a given statistic changed between two cases
% (typically "before" and "after" measurements).
%
% NOTE - This kludges a box chart into showing -2sigma/mean/2sigma rather
% than quartiles and median.
%
% FIXME - This expects 2d data, not 3d. The caller will have to select the
% slice corresponding to the desired trial.
%
% "meanbefore" is a matrix indexed by (firstidx,secondidx) with mean values.
% "confbefore" is a matrix indexed by (firstidx,secondidx) with confidence
%   interval radii (positive values).
% "meanafter" is a matrix indexed by (firstidx,secondidx) with mean values.
% "confafter" is a matrix indexed by (firstidx,secondidx) with confidence
%   interval radii (positive values).
% "firstlabels" is a cell array with the first series of channel labels.
% "secondlabels" is a cell array with the second series of channel labels.
% "yrange" [ min max ] is the axis range for statistics values.
% "ycursors" is a vector containing y values to render horizontal cursors on.
% "ytitle" is the label applied to the statistics axis.
% "figtitle" is the title to use for the figure.
% "fname" is the filename to use when saving the figure.


% Get metadata.

firstcount = size(meanbefore,1);
secondcount = size(meanbefore,2);

pairmask = nlUtil_getPairMask(firstlabels, secondlabels);



%
% Prepare box chart data.

% FIXME - Kludge. Build fake data series for box plots.
% Using [ (u-conf) (u-conf) u (u+conf) (u+conf) ] gets confidence boxes.

% "bin" is the channel pair.
% "set" is before/after.

datavalues = [];
databinvalues = {};
datasetlabels = {};

for firstidx = 1:firstcount
  for secondidx = 1:secondcount
    if pairmask(firstidx,secondidx)

      % Pair label.
      thislabel = [ firstlabels{firstidx} ' - ' secondlabels{secondidx} ];

      % "Before" data.

      thismean = meanbefore(firstidx,secondidx);
      thisconf = confbefore(firstidx,secondidx);
      thislower = thismean - thisconf;
      thisupper = thismean + thisconf;

      thisdata = [ thislower thislower thismean thisupper thisupper ];
      thisbins = cell(size(thisdata));
      thisbins(:) = { thislabel };
      thissets = cell(size(thisdata));
      thissets(:) = { 'before' };

      datavalues = [ datavalues thisdata ];
      databinvalues = [ databinvalues thisbins ];
      datasetlabels = [ datasetlabels thissets ];

      % "After" data.

      thismean = meanafter(firstidx,secondidx);
      thisconf = confafter(firstidx,secondidx);
      thislower = thismean - thisconf;
      thisupper = thismean + thisconf;

      thisdata = [ thislower thislower thismean thisupper thisupper ];
      thisbins = cell(size(thisdata));
      thisbins(:) = { thislabel };
      thissets = cell(size(thisdata));
      thissets(:) = { 'after' };

      datavalues = [ datavalues thisdata ];
      databinvalues = [ databinvalues thisbins ];
      datasetlabels = [ datasetlabels thissets ];

    end
  end
end



%
% Render the box charts.

cols = nlPlot_getColorPalette();
ycursordefs = {};
for cidx = 1:length(ycursors)
  ycursordefs{cidx} = { ycursors(cidx), cols.gry, '' };
end

euPlot_plotMultipleBoxCharts( datavalues, databinvalues, datasetlabels, ...
  'linear', 'linear', [], yrange, '', ytitle, false, 0.5, ...
  {}, ycursordefs, 'northeastoutside', figtitle, fname );



% Done.
end


%
% This is the end of the file.
