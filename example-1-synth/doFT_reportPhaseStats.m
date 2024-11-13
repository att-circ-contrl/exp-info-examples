function reporttext = doFT_reportPhaseStats( phasedata )

% function reporttext = doFT_reportPhaseStats( phasedata )
%
% This compiles a human-readable report of a statistical analysis of
% phase alignment between pairs of channels across trials and time windows.
%
% "phasedata" is a structure returned by eiCalc_calcTrialPhaseStats().
%
% "reporttext" is a human-readable summary of relevant statistics.


reporttext = '';


% FIXME - Just reporting aggregate stats, not pairwise.


%
% Phase difference.

% Convert to degrees.
thisdata = phasedata.phasedifftrialsdata * 180 / pi;

[ thismean thisdev scratch ] = ...
  nlProc_calcMatrixStats( thisdata, [ 10 25 50 75 90 ] );

reporttext = [ reporttext sprintf( ...
  'Angle:    mean %d  sd %d      median %d  iqr %d  idr %d\n', ...
  round(thismean), round(thisdev), round(scratch(3)), ...
  round(scratch(4) - scratch(2)), round(scratch(5) - scratch(1)) ) ];


%
% Phase lock value.

thisdata = phasedata.plvtrialsdata;

[ thismean thisdev scratch ] = ...
  nlProc_calcMatrixStats( thisdata, [ 10 25 50 75 90 ] );

reporttext = [ reporttext sprintf( ...
  '  PLV:    mean %.2f  sd %.2f      median %.2f  iqr %.2f  idr %.2f\n', ...
  thismean, thisdev, scratch(3), ...
  scratch(4) - scratch(2), scratch(5) - scratch(1) ) ];



% Done.
end


%
% This is the end of the file.
