% Common setup tasks.
% Written by Christopher Thomas.
%
% This sets up paths, initializes Open Ephys, turns off warnings, etc.
% I'm tired of constantly having to retype this.


% Paths.

addpath('../libraries/lib-entropy');
addpath('../libraries/lib-exp-info');
addpath('../libraries/lib-exp-utils');
addpath('../libraries/lib-looputil');
addpath('../libraries/lib-fieldtrip');

%addpath('../libraries/lib-openephys');
%addpath('../libraries/lib-npy-matlab');

addPathsExpUtils;
addPathsExpInfo;
addPathsLoopUtil;


% Matlab warnings.

oldwarnstate = warning('off');


% Field Trip initialization and messages.

evalc('ft_defaults');

ft_notice('off');
ft_info('off');
ft_warning('off');


% LoopUtil configuration.

% We use about 1 GB per channel-hour. The sweet spot to reduce loading time
% is 4 channels or more.
nlFT_setMemChans(4);


%
% This is the end of the file.
