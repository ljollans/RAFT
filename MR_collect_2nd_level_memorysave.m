function MR_collect_2nd_level_memorysave(design)
% collects and combines intermittent files saved in RAFT_2nd_level
%
% for comments and questions please contact lee.jollans@gmail.com

cd(design.saveto);
design.nvars=size(design.data,2);
fprintf('Aggregating Middle Fold Analysis and deleting intermediate files\n', size(design.data));
%%
for mainfold=1:design.numFolds
    eval(['LHmerit' num2str(mainfold) '=NaN(design.numFolds, design.numFolds, size(design.merit_thresholds,2));']);
    eval(['vars2use' num2str(mainfold) '=zeros(design.numFolds, design.numFolds, size(design.merit_thresholds,2),size(design.data,2));']);
    for subfolds=1:design.numFolds
        load([design.saveto filesep 'Second_level_results_' num2str(mainfold) '_' num2str(subfolds)]);
        try
        eval(['LHmerit' num2str(mainfold) '(subfolds,:,:)=(LHmerit2save);']);
        eval(['vars2use' num2str(mainfold) '(subfolds,:,:,:)=vars2use2save;']);
        catch
            pause
        end
        end
    eval(['vars2use' num2str(mainfold) '=logical(vars2use' num2str(mainfold) ');']);
    fname=['LHmerit' num2str(mainfold)];  save(fname, fname);
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