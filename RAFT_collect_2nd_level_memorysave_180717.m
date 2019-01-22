function RAFT_collect_2nd_level_memorysave(design)
% collects and combines intermittent files saved in RAFT_2nd_level
%
% for comments and questions please contact lee.jollans@gmail.com

cd(design.saveto);
design.nvars=size(design.data,2);
fprintf('Aggregating Middle Fold Analysis and deleting intermediate files\n', size(design.data));
%%
for mainfold=1:design.numFolds
switch(design.type)
case 'linear'
eval(['LHmerit' num2str(mainfold) '=NaN(design.numFolds, design.numFolds, size(design.merit_thresholds,2),design.numAlphas, design.numLambdas);']);
case 'logistic'
eval(['LHmerit' num2str(mainfold) '.accuracy=NaN(design.numFolds, design.numFolds, size(design.merit_thresholds,2),design.numAlphas, design.numLambdas);']);
eval(['LHmerit' num2str(mainfold) '.sensitivity=NaN(design.numFolds, design.numFolds, size(design.merit_thresholds,2),design.numAlphas, design.numLambdas);']);
eval(['LHmerit' num2str(mainfold) '.specificity=NaN(design.numFolds, design.numFolds, size(design.merit_thresholds,2),design.numAlphas, design.numLambdas);']);
eval(['LHmerit' num2str(mainfold) '.overall_AUC=NaN(design.numFolds, design.numFolds, size(design.merit_thresholds,2),design.numAlphas, design.numLambdas);']);
eval(['LHmerit' num2str(mainfold) '.F1score=NaN(design.numFolds, design.numFolds, size(design.merit_thresholds,2),design.numAlphas, design.numLambdas);']);
end
    eval(['vars2use' num2str(mainfold) '=zeros(design.numFolds, design.numFolds, size(design.merit_thresholds,2),design.numAlphas, design.numLambdas,size(design.data,2));']);
    eval(['lambda_values' num2str(mainfold) '=NaN(design.numFolds, design.numFolds, size(design.merit_thresholds,2),design.numAlphas, design.numLambdas);']);
    for subfolds=1:design.numFolds
        load([design.saveto filesep 'Second_level_results_' num2str(mainfold) '_' num2str(subfolds)]);
switch(design.type)
case 'linear',
        eval(['LHmerit' num2str(mainfold) '(subfolds,:,:,:,:)=LHmerit2save;']);
case 'logistic'
eval(['LHmerit' num2str(mainfold) '.accuracy(subfolds,:,:,:,:)=LHmerit2save.accuracy;']);
eval(['LHmerit' num2str(mainfold) '.sensitivity(subfolds,:,:,:,:)=LHmerit2save.sensitivity;']);
eval(['LHmerit' num2str(mainfold) '.specificity(subfolds,:,:,:,:)=LHmerit2save.specificity;']);
eval(['LHmerit' num2str(mainfold) '.overall_AUC(subfolds,:,:,:,:)=LHmerit2save.overall_AUC;']);
eval(['LHmerit' num2str(mainfold) '.F1score(subfolds,:,:,:,:)=LHmerit2save.F1score;']);
end
        eval(['vars2use' num2str(mainfold) '(subfolds,:,:,:,:,:)=vars2use2save;']);
        eval(['lambda_values' num2str(mainfold) '(subfolds,:,:,:,:)=lambda_values2save;']);
    end
    eval(['vars2use' num2str(mainfold) '=logical(vars2use' num2str(mainfold) ');']);
    fname=['LHmerit' num2str(mainfold)];  save(fname, fname);
    fname=['lambda_values' num2str(mainfold)];  save(fname, fname);
    fname=['vars2use' num2str(mainfold)];
    try
      save(fname, fname);
    catch
      save(fname, fname, '-v7.3');
    end
end

for folds=1:(design.numFolds*design.numFolds)
        [mainfold subfolds]=ind2sub([design.numFolds design.numFolds], folds);  
        delete([design.saveto filesep 'Second_level_results_' num2str(mainfold) '_' num2str(subfolds) '.mat']);
end
end
