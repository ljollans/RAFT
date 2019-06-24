# regularized adaptive feature thresholding (RAFT)
    Copyright (C) 2019,  Lee Jollans

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.


[![DOI](https://zenodo.org/badge/74484393.svg)](https://zenodo.org/badge/latestdoi/74484393)

Basic use:

create design files using create_design.m

run a model using all_ML.m

randomeffects2.m was used to create simulated data files

simlooper.m is an example of how multiple different model types might be run while simulating data on the fly

looper.m is a function to run multiple iterations of the same model

wdhlooper_aggregate.m summarizes resutls from multiple model iterations

Some of the code used here comes from other sources. Authors will always be listed in the .m file itself.
To run the majority of this code you will also need to download the following:

% boostrapal.m (included in this repository)

% glmnet - http://web.stanford.edu/~hastie/glmnet_matlab/
(note that use of glmnet in MATLAB on MacOS does not always work. An alternative implementation of glmnet is available in R. A work-around including a wrapper script in MATLAB incorporating glmnet in R has been implemented by ljollans. Get in touch for more details.)

% fastAUC - http://www.mathworks.com/matlabcentral/fileexchange/42860-fast-auc-calculator-and-roc-curve-plotter/content/fastAUC.m


NOTE FOR IOS USERS: glmnet does not work in many macos and linux instances of MATLAB at the moment. We have implemented a simple wrapper for the elastic net in R (multinomial) and python (linear and multinomial). Upload of these wrappers to follow. Further implementation of result interpretation etc in python in progress.
