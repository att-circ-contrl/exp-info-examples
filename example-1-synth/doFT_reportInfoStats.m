function reporttext = doFT_reportInfoStats( ...
  ampbefore, ampafter, lagbefore, lagafter, ampsigmask, lagsigmask, ...
  destlabels, srclabels )

% function reporttext = doFT_reportInfoStats( ...
%   ampbefore, ampafter, lagbefore, lagafter, ampsigmask, lagsigmask, ...
%   destlabels, srclabels )
%
% This compiles a human-readable report of statistically significant
% changes in before-and-after pairwise similarity peak amplitude and time lag.
%
% Reporting is suppressed for cells with NaN values.
%
% "ampbefore" is a matrix indexed by (destidx,srcidx,trialidx) with mean
%   pairwise similarity amplitude values from the "before" period.
% "ampafter" is a matrix indexed by (destidx,srcidx,trialidx) with mean
%   pairwise similarity amplitude values from the "after" period.
% "lagbefore" is a matrix indexed by (destidx,srcidx,trialidx) with mean
%   pairwise similarity lag values from the "before" period.
% "lagafter" is a matrix indexed by (destidx,srcidx,trialidx) with mean
%   pairwise similarity lag values from the "after" period.
% "ampsigmask" is a matrix indexed by (destidx,srcidx,trialidx) that's true
%   for significant amplitude changes and false otherwise.
% "lagsigmask" is a matrix indexed by (destidx,srcidx,trialidx) that's true
%   for significant time lag changes and false otherwise.
% "destlabels" is a cell array with the destination channel labels.
% "srclabels" is a cell array with the source channel labels.
%
% "reporttext" is a character vector with a human-readable summary of the
%   significant changes between "before" and "after" measurements.


% Get metadata.

destcount = length(destlabels);
srccount = length(srclabels);

trialcount = size(ampbefore,3);


% Mask off self-comparison, permutations, and NaN data.

pairmask = nlUtil_getPairMask(destlabels, srclabels);

% Augment this to match data dimensions.
for tidx = 2:trialcount
  pairmask(:,:,tidx) = pairmask(:,:,1);
end

goodmask = isnan(ampbefore) | isnan(ampafter) ...
  | isnan(lagbefore) | isnan(lagafter);


% FIXME - Debugging. Show NaN pairs.
goodmask = pairmask;


%
% Build a human-readable report.

reporttext = '';

newline = sprintf('\n');

for destidx = 1:destcount
  for srcidx = 1:srccount
    if any( goodmask(destidx,srcidx,:) )

      thisline = [ destlabels{destidx} ' from ' srclabels{srcidx} ':' ];

      % FIXME - This will get spammy if we have lots of trials!
      for tidx = 1:trialcount
        if goodmask(destidx,srcidx,tidx)

          % FIXME - We don't have the real trial numbers!
          thisline = [ thisline sprintf('\n  %Tr%04d:    ', tidx) ];

          thisfragment = '';
          if ampsigmask(destidx, srcidx)
            thisfragment = sprintf( 'amp %+.2f to %+.2f', ...
              ampbefore(destidx,srcidx,tidx), ...
              ampafter(destidx,srcidx,tidx) );
          else
            % FIXME - No way to determine which should get priority.
            thisamp = mean( [ ampbefore(destidx,srcidx,tidx), ...
              ampafter(destidx,srcidx,tidx) ] );

            thisfragment = sprintf( '(amp ns) %+.2f', thisamp );
          end
          thisline = [ thisline sprintf( '%-18s    ', thisfragment ) ];

          if lagsigmask(destidx, srcidx)
            thisline = [ thisline sprintf( 'lag %+d to %+d', ...
              round( lagbefore(destidx,srcidx,tidx) ), ...
              round( lagafter(destidx,srcidx,tidx) ) ) ];
          else
            % FIXME - No way to determine which should get priority.
            thislag = mean( [ lagbefore(destidx,srcidx,tidx), ...
              lagafter(destidx,srcidx,tidx) ] );

            thisline = [ thisline sprintf( '(lag ns) %+d', round(thislag) ) ];
          end

          reporttext = [ reporttext thisline newline ];

        end
      end

    end
  end
end



% Done.
end


%
% This is the end of the file.
