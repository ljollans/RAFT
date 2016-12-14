function [LHmerit, vars2use, lambda_values] = RAFT_collect_2nd_level(design)
% collects and combines intermittent files saved in RAFT_2nd_level
%
% for comments and questions please contact lee.jollans@gmail.com

cd(design.saveto);

fprintf('Aggregating Middle Fold Analysis and deleting intermediate files\n', size(design.data));
%%
LHmerit=NaN(design.numFolds, design.numFolds, design.numFolds, size(design.merit_thresholds,2),design.numAlphas, design.numLambdas);
vars2use=zeros(design.numFolds, design.numFolds, design.numFolds, size(design.merit_thresholds,2),design.numAlphas, design.numLambdas,design.nvars);
lambda_values=NaN(design.numFolds, design.numFolds, design.numFolds, size(design.merit_thresholds,2),design.numAlphas, design.numLambdas);

for folds=1:(design.numFolds*design.numFolds)
        [mainfold subfolds]=ind2sub([design.numFolds design.numFolds], folds);  
        load([design.saveto filesep 'Second_level_results_' num2str(mainfold) '_' num2str(subfolds)]);
        LHmerit(mainfold, subfolds,:,:,:,:)=LHmerit2save;
        vars2use(mainfold,subfolds,:,:,:,:,:)=vars2use2save;
        lambda_values(mainfold,subfolds,:,:,:,:)=lambda_values2save;
end

vars2use=logical(vars2use);

save('LHmerit', 'LHmerit')
try
save('vars2use', 'vars2use')
catch
save('vars2use', 'vars2use', '-v7.3')
end
save('lambda_values', 'lambda_values')
for folds=1:(design.numFolds*design.numFolds)
        [mainfold subfolds]=ind2sub([design.numFolds design.numFolds], folds);  
        delete([design.saveto filesep 'Second_level_results_' num2str(mainfold) '_' num2str(subfolds) '.mat']);
end
end