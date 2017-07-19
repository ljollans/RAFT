function MR_collect_2nd_level(design)
% collects and combines intermittent files saved in RAFT_2nd_level
%
% for comments and questions please contact lee.jollans@gmail.com

cd(design.saveto);
design.nvars=size(design.data,2);
fprintf('Aggregating Middle Fold Analysis and deleting intermediate files\n', size(design.data));
%%
LHmerit=NaN(design.numFolds, design.numFolds, design.numFolds, size(design.merit_thresholds,2));
vars2use=zeros(design.numFolds, design.numFolds, design.numFolds, size(design.merit_thresholds,2),size(design.data,2));

for folds=1:(design.numFolds*design.numFolds)
        [mainfold subfolds]=ind2sub([design.numFolds design.numFolds], folds);  
        load([design.saveto filesep 'Second_level_results_' num2str(mainfold) '_' num2str(subfolds)]);
        if size(LHmerit2save,1)==100
            LHmerit(mainfold, subfolds,:,:)=LHmerit2save{folds};
        else
            LHmerit(mainfold, subfolds,:,:)=LHmerit2save;
        end
        if size(vars2use2save,1)==100
            vars2use(mainfold,subfolds,:,:,:)=vars2use2save{folds};
        else
            vars2use(mainfold,subfolds,:,:,:)=vars2use2save;
        end        
end

vars2use=logical(vars2use);

save('LHmerit', 'LHmerit')
try
save('vars2use', 'vars2use')
catch
save('vars2use', 'vars2use', '-v7.3')
end

for folds=1:(design.numFolds*design.numFolds)
        [mainfold subfolds]=ind2sub([design.numFolds design.numFolds], folds);  
        delete([design.saveto filesep 'Second_level_results_' num2str(mainfold) '_' num2str(subfolds) '.mat']);
end
end