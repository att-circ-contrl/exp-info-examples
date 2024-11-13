function doFT_plotPairWaves( ...
  pairdata, firstlabels, secondlabels, trialnums, ...
  xseries, smoothsamples, yrange, ycursors, xtitle, ytitle, figtitle, fbase )

% function doFT_plotPairWaves( ...
%   pairdata, firstlabels, secondlabels, trialnums, ...
%   xseries, smoothsamples, yrange, ycursors, xtitle, ytitle, figtitle, fbase )
%
% This plots several pairwise data curves.
%
% "pairdata" is a matrix indexed by (firstidx, secondidx, trialidx, xidx).
% "firstlabels" is a cell array with the first series of channel labels.
% "secondlabels" is a cell array with the second series of channel labels.
% "trialnums" is a vector containing trial numbers (arbitrary integer values
%   used as labels for each trial).
% "xseries" is a vector with independent variable values.
% "smoothsamples" is the window size for smoothing data in the x direction,
%   or 0 or NaN to not smooth.
% "yrange" [ min max ] is the plot range for data values, or [] to auto-range.
% "ycursors" is a vector containing Y values to draw horizontal cursors at.
% "xtitle" is the label applied to the X axis.
% "ytitle" is the label applied to the Y axis.
% "figtitle" is the title to use for the figure.
% "fbase" is the prefix to use when building plot filenames.
%
% No return value.


firstcount = size(pairdata,1);
secondcount = size(pairdata,2);
trialcount = size(pairdata,3);
xcount = size(pairdata,4);


want_smooth = true;
if isnan(smoothsamples) || (smoothsamples < 1)
  want_smooth = false;
end



% Get the number of pairs.
% Also get a mask matrix. This is intended to work for partly-overlapping
% channel lists.

pairmask = nlUtil_getPairMask(firstlabels, secondlabels);
paircount = sum(sum(pairmask));



% Get a palette.
% FIXME - Hardcoding colours.

cols = nlPlot_getColorPalette();
paircols = nlPlot_getColorSpread( cols.brn, paircount, 260 );



% Set up the figure and render the plots.

thisfig = figure();
figure(thisfig);


for trialidx = 1:trialcount

  trialtitle = '';
  triallabel = '';

  if trialcount > 1
    trialtitle = sprintf( ' - Tr %04d', trialnums(trialidx) );
    triallabel = sprintf( '-tr%04d', trialnums(trialidx) );
  end


  % Render this plot's curves.

  clf('reset');
  hold on;

  pairidx = 0;
  for firstidx = 1:firstcount
    for secondidx = 1:secondcount
      if pairmask(firstidx,secondidx)

        pairidx = pairidx + 1;
        thiscol = paircols{pairidx};

        thisdata = pairdata(firstidx, secondidx, trialidx, :);
        thisdata = reshape( thisdata, size(xseries) );

        if want_smooth
          % Triangular window with about 1.5x the specified size.
          thisdata = movmean(thisdata, smoothsamples);
          thisdata = movmean(thisdata, smoothsamples);
        end

        thislabel = [ firstlabels{firstidx} ' - ' secondlabels{secondidx} ];

        plot( xseries, thisdata, 'Color', thiscol, 'DisplayName', thislabel );

      end
    end
  end

  % Add cursors.
  xmin = min(xseries);
  xmax = max(xseries);
  for cidx = 1:length(ycursors)
    thisy = ycursors(cidx);
    plot( [ xmin xmax ], [ thisy thisy ], 'Color', cols.gry, ...
      'HandleVisibility', 'off' );
  end

  hold off;

  if ~isempty(yrange)
    ylim(yrange);
  end

  xlabel(xtitle);
  ylabel(ytitle);

  title(figtitle, trialtitle);

  legend('Location', 'northeastoutside');

  saveas(thisfig, [ fbase triallabel '.png' ]);

end


% Dispose of the figure.
close(thisfig);



% Done.
end


%
% This is the end of the file.
