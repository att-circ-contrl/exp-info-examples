# Connectivity Analysis Example Scripts

## Overview

This is a set of example scripts that use the `exp-info` library to analyze
connectivity and information transfer between channels in neural recordings.

These are not well-polished tutorial scripts. These are derived from Chris's
testing and pilot analysis code. They should still be useful as a starting
point for new analysis scripts.


## Libraries

The scripts need libraries from the following ACC Lab GitHub repositories:

* `ConditionalEntropy` (low-level information calculations):
<https://github.com/att-circ-contrl/ConditionalEntropy>
* `exp-info` (high-level information calculations)
<https://github.com/att-circ-contrl/exp-info>
* `exp-preproc` (pre-processing and iteration routines)
<https://github.com/att-circ-contrl/exp-preproc>
* `exp-utils` (miscellaneous support functions):
<https://github.com/att-circ-contrl/exp-utils>
* `LoopUtil` (miscellaneous support functions):
<https://github.com/att-circ-contrl/LoopUtil>

The following external libraries are also needed:

* The Field Trip toolbox:
<https://github.com/fieldtrip/fieldtrip>
* The Open Ephys "`analysis tools`" library (for reading data saved in Open
Ephys's format):
<https://github.com/open-ephys/analysis-tools>
* The `npy-matlab` library (needed by Open Ephys's library; it comes with it).


## Folders

* `datasets` --
Abbreviated datasets read by the example scripts.
* `example-1-synth` --
Analysis script examining a synthetic dataset with known communication
patterns (generated using the `SynthRobinson` repository). This looks for
"routing states" (pairs of channels with directional information flow).
* `libraries` --
Symbolic links to the various libraries used by the scripts. These will only
work on my system; modify them or add library paths to Matlab manually to
run the scripts on your own system.



*(This is the end of the file.)*
