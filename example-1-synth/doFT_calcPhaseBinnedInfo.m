function [ vsphasetime vsphaselagboth vsphaselagafter ] = ...
  doFT_calcPhaseBinnedInfo( ftdata_dest, ftdata_src, ...
    win_params, infofuncphase, datafield, minplv, ...
    phasebincount, phasebinwidth_deg, timesmooth_ms, afterspan_ms )

% function [ vsphasetime vsphaselagboth vsphaselagafter ] = ...
%   doFT_calcPhaseBinnedInfo( ftdata_dest, ftdata_src, ...
%     win_params, infofuncphase, datafield, minplv, ...
%     phasebincount, phasebinwidth_deg, timesmooth_ms, afterspan_ms )
%
% This calculates average information similarity vs phase bin and lag and vs
% phase bin and time. Phase is the average phase difference between the
% two input signals during the time window being evaluated.
%
% "ftdata_dest" is a ft_datatype_raw structure with trial data for the
%   putative destination channels.
% "ftdata_src" is a ft_datatype_raw structure with trial data for the
%   putative source channels.
% "win_params" is a structure giving time window and time lag information,
%   per TIMEWINLAGSPEC.txt.
% "infofuncphase" is a function handle with the form:
%   data = infofuncphase( ftdata_dest, ftdata_src, winparams, phaseparams )
% "datafield" is the base name of the data field produced by "infofuncphase".
%   This is used as a prefix for building field names, per TIMEWINLAGDATA.txt.
% "minplv" is the minimum phase-lock value required for phase binning.
%   Phase-lock value is (1 - circular variance).
% "phasebincount" is the number of phase bins to test.
% "phasebinwidth_deg" is the width of the phase bins, in degrees, or NaN to
%   use 360/count.
% "timesmooth_ms" is the width of the smoothing window to use for window
%   times, in milliseconds, or NaN to not smooth.
% "afterspan_ms" [ min max ] is the time range to use for "after simulation"
%   time-averaging.
%
% "vsphasetime" is a structure with the following fields:
%   "destchans" is a cell array with FT channel names for the putative
%     destination channels.
%   "srcchans" is a cell array with FT channel names for the putative
%     source channels.
%   "phaselist_deg" is a vector containing phase bin midpoints in degrees.
%   "windowlist_ms" is a vector containing timestamps in milliseconds
%     specifying where the middle of each analysis window is.
%   "trialnums" is a vector containing trial numbers (arbitrary integer
%     values used as labels for each trial).
%   "avg" is a matrix indexed by (destchan, srcchan, trial, winidx, phase)
%     containing information similarity values averaged across time lag.
%   "dev" is a matrix indexed by (destchan, srcchan, trial, winidx, phase)
%     containing the standard deviation of the information similarity averages.
%
% "vsphaselagboth" is a structure with the following fields:
%   "destchans" is a cell array with FT channel names for the putative
%     destination channels.
%   "srcchans" is a cell array with FT channel names for the putative
%     source channels.
%   "phaselist_deg" is a vector containing phase bin midpoints in degrees.
%   "delaylist_ms" is a vector containing the time lags tested in
%     milliseconds.
%   "trialnums" is a vector containing trial numbers (arbitrary integer
%     values used as labels for each trial).
%   "avg" is a matrix indexed by (destchan, srcchan, trial, lagidx, phase)
%     containing information similarity values averaged across window time.
%   "dev" is a matrix indexed by (destchan, srcchan, trial, lagidx, phase)
%     containing the standard deviation of the information similarity averages.
%
% "vsphaselagafter" is a structure with the same fields as "vsphaselagboth",
%   but only averaging across time windows within the "after" time span.


%
% Get phase bin locations.

phaselist_deg = 1:phasebincount;

% Put zero in the middle of the list.
phaselist_deg = phaselist_deg - round( (0.5 * phasebincount) + 0.1 );

phaselist_deg = 360 * (phaselist_deg / phasebincount);


%
% Walk through phase bins, building output.

vsphasetime = struct([]);
vsphaselag = struct([]);

phaseparams = struct();
phaseparams.phasetargetdeg = 0;
phaseparams.acceptwidthdeg = phasebinwidth_deg;
phaseparams.minplv = minplv;

for pidx = 1:phasebincount

  phaseparams.phasetargetdeg = phaselist_deg(pidx);

  % FIXME - Diagnostics.
  timescratch = tic;

  timelagdata = ...
    infofuncphase( ftdata_dest, ftdata_src, win_params, phaseparams );

  % FIXME - Diagnostics.
  durstring = nlUtil_makePrettyTime(toc(timescratch));
  disp(sprintf( '.. Calculated for %d degrees in %s.', ...
    round( phaselist_deg(pidx) ), durstring ));

  % This recognizes the NaN "don't smooth" case, so no special-casing needed.
%  smooth_type = 'coarse';
  smooth_type = 'smooth';
  timelagdata = eiCalc_smoothTimeLagAverages( ...
    timelagdata, { datafield }, timesmooth_ms, NaN, smooth_type );

  [ thisvstime thisvslagafter ] = eiCalc_collapseTimeLagAverages( ...
    timelagdata, datafield, {afterspan_ms}, {[]} );
  [ thisvstime thisvslagboth ] = eiCalc_collapseTimeLagAverages( ...
    timelagdata, datafield, {[]}, {[]} );

  if isempty(vsphasetime)

    % Initialize output.

    vsphasetime = thisvstime;
    vsphasetime = rmfield(vsphasetime, 'delayrange_ms');
    vsphasetime.('phaselist_deg') = phaselist_deg;

    vsphaselagboth = thisvslagboth;
    vsphaselagboth = rmfield(vsphaselagboth, 'windowrange_ms');
    vsphaselagboth.('phaselist_deg') = phaselist_deg;

    vsphaselagafter = thisvslagafter;
    vsphaselagafter = rmfield(vsphaselagafter, 'windowrange_ms');
    vsphaselagafter.('phaselist_deg') = phaselist_deg;

  else

    % Append output.

    vsphasetime.avg(:,:,:,:,pidx) = thisvstime.avg;
    vsphasetime.dev(:,:,:,:,pidx) = thisvstime.dev;

    vsphaselagboth.avg(:,:,:,:,pidx) = thisvslagboth.avg;
    vsphaselagboth.dev(:,:,:,:,pidx) = thisvslagboth.dev;

    vsphaselagafter.avg(:,:,:,:,pidx) = thisvslagafter.avg;
    vsphaselagafter.dev(:,:,:,:,pidx) = thisvslagafter.dev;

  end

end


% Done.
end

%
% This is the end of the file.
