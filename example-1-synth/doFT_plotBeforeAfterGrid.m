function doFT_plotBeforeAfterGrid( ...
  ampbefore, ampafter, lagbefore, lagafter, ...
  flagsigamp, flagsiglag, firstlabels, secondlabels, ...
  colmapamp, colmaplag, infotitle, infoabbrv, titleprefix, fname )

% function doFT_plotBeforeAfterGrid( ...
%   ampbefore, ampafter, lagbefore, lagafter, ...
%   flagsigamp, flagsiglag, firstlabels, secondlabels, ...
%   colmapamp, colmaplag, infotitle, infoabbrv, titleprefix, fname )
%
% This makes a multi-pane plot showing significant changes in the amplitude
% and time-lag of the pairwise information peak between two cases (typically
% "before" and "after" measurements).
%
% FIXME - This expects 2d data, not 3d. The caller will have to select the
% slice corresponding to a desired trial.
%
% This will squash any cells that are self-comparison or permutations of
% pairs already plotted.
% This assumes that invalid cells have already been NaN-squashed.
%
% "ampbefore" is a matrix indexed by (firstidx,secondidx) with
%    pairwise information peak amplitude values.
% "ampafter" is a matrix indexed by (firstidx,secondidx) with
%    pairwise information peak amplitude values.
% "lagbefore" is a matrix indexed by (firstidx,secondidx) with lag values.
% "lagafter" is a matrix indexed by (firstidx,secondidx) with lag values.
% "flagsigamp" is a matrix indexed by (firstidx,secondidx) that's true for
%   cells that have significant changes in amplitude and false otherwise.
% "flagsiglag" is a matrix indexed by (firstidx,secondidx) that's true for
%   cells that have significant changes in time lag and false otherwise.
% "firstlabels" is a cell array with the first series of channel labels.
% "secondlabels" is a cell array with the second series of channel labels.
% "colmapamp" is a colour map used for plotting pairwise information amplitude.
% "colmaplag" is a colour map used for plotting pairwise information time lag.
% "infotitle" is the name of the data value being plotted (e.g.
%   'Cross-Correlation').
% "infoabbrv" is an abbreviation of the name of the data value (e.g. 'XC').
% "titleprefix" is a character vector used for building plot and subplot
%   titles.
% "fname" is the name of the file to write to.
%
% No return value.


% Get metadata.

firstcount = length(firstlabels);
secondcount = length(secondlabels);


% Squash self-comparisons and permutations.

pairmask = nlUtil_getPairMask(firstlabels, secondlabels);

ampbefore(~pairmask) = NaN;
ampafter(~pairmask) = NaN;
lagbefore(~pairmask) = NaN;
lagafter(~pairmask) = NaN;

flagsigamp(~pairmask) = false;
flagsiglag(~pairmask) = false;


% Get data ranges.

if true
  % Range to cover the data.

  scratch = max( abs(ampbefore), abs(ampafter) );
  scratch = max( scratch, [], 'all' );
  scratch = max(scratch, 1e-30);
  amprange = [ (-scratch) scratch ];

  scratch = max( abs(lagbefore), abs(lagafter) );
  scratch = max( scratch, [], 'all' );
  scratch = max(scratch, 1e-30);
  lagrange = [ (-scratch) scratch ];
else
  % FIXME - Hard-coded ranges for consistency.
  amprange = [ -1.0 1.0 ];
  lagrange = [ -20 20 ];
end



%
% Set up plotting.

thisfig = figure();

figure(thisfig);
clf('reset');

% Adjust geometry so that subplots are legible.

oldpos = thisfig.Position;
oldwidth = oldpos(3);
oldheight = oldpos(4);

% A rectangular tile fits the colourbar better.
tilewidth = max(oldwidth, oldheight);
tileheight = min(oldwidth, oldheight);

% We're making a 2x2 grid of plots.
plotrows = 2;
plotcols = 2;

newpos = oldpos;
newpos(3) = tilewidth * plotcols;
newpos(4) = tileheight * plotrows;

thisfig.Position = newpos;



%
% Render the plots.

plotgrid = { ampbefore, ampafter ; lagbefore, lagafter };
sigbyrow = { flagsigamp ; flagsiglag };
mapbyrow = { colmapamp ; colmaplag };
rangebyrow = { amprange ; lagrange };

rowlabels = { [ infoabbrv 'Peak Amp' ] ; [ infoabbrv ' Time Lag' ] };
collabels = { 'Before', 'After' };
rowscales = { infotitle ; 'Time Lag (ms)' };

for plotrowidx = 1:plotrows
  for plotcolidx = 1:plotcols

    % Get this pane's data and decoration info.

    thisplotdata = plotgrid{plotrowidx,plotcolidx};
    thisplotsig = sigbyrow{plotrowidx};

    thiscolmap = mapbyrow{plotrowidx};
    thiscolscale = rowscales{plotrowidx};
    thisrange = rangebyrow{plotrowidx};

    thistitle = [ rowlabels{plotrowidx} ' ' collabels{plotcolidx} ];
    % NOTE - If we're using a global title, we don't need to prepend this.
%    thistitle = [ titleprefix ' - ' thistitle ];


    % Select this pane.

    tileidx = (plotrowidx - 1) * plotcols;
    tileidx = tileidx + (plotcolidx - 1);
    tileidx = tileidx + 1;

    subplot( plotrows, plotcols, tileidx );
    thisax = gca;


    % Render this pane.

    euPlot_axesPlotSignificantHeatmap( thisax, thisplotdata, thisplotsig, ...
      secondlabels, firstlabels, [], [], 'Destination', 'Source', thistitle );


    % Flip to match raster order for rendering channels.
    set(thisax, 'YDir', 'reverse');

    colormap(thisax, thiscolmap);
    clim(thisax, thisrange);

    thiscol = colorbar(thisax);
    thiscol.Label.String = thiscolscale;

  end
end

% Add the global title.

titletext = sgtitle(thisfig, titleprefix);
titletext.FontSize = 20;


% Finished rendering.

saveas(thisfig, fname);



%
% Finished with this figure.

thisfig.Position = oldpos;
close(thisfig);


% Done.
end


%
% This is the end of the file.
