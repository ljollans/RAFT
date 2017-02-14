function [params2pick, Merit, Beta, Vars2pick_main, stats]=RAFT(design)
% Regularized Adaptive Feature Thresholding (RAFT)
% 3-tier Machine Learning, Elastic Net with Feature Selection
% Authors: Lee Jollans, Robert Whelan, Richard Watts
%
% create_design.m - this may help in constructing the design file
%% Files used in this code:
% boostrapal.m - Creates datasets for bagging
% glmnet - http://web.stanford.edu/~hastie/glmnet_matlab/
% fastAUC - http://www.mathworks.com/matlabcentral/fileexchange/42860-fast-auc-calculator-and-roc-curve-plotter/content/fastAUC.m

%%
cd(design.saveto)
design.nvars=size(design.data,2);
[merit_per_var] = RAFT_FS(design); 
[design, pass_vars]=RAFT_do_thresh(design, merit_per_var);
RAFT_2nd_Level(design, pass_vars, [1:(design.numFolds*design.numFolds)]);
if size(design.data,2)>=1000
    RAFT_collect_2nd_level_memorysave(design);
else
    RAFT_collect_2nd_level(design);
end
[params2pick, Merit, Beta, Vars2pick_main,  GetProbs, design, stats] = RAFT_Model_Selection_140217(design, 1);
end