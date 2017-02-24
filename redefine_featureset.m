function redefine_featureset(olddir, newdir, numreps, vars2keep, deleteprev)
% function to re-run an analysis with a subset of the features which were
% used in a previous analysis without having to redo teh feature thresholding step
% example usage: redefine_featureset('/home/lee/local/MS ML/Results_CMSMLlinCOGNzMO12_P3bVis',...
% '/home/lee/local/MS ML/Results_CMSMLlinCOGNzMO12_P3bVis_new', 20, 2:1620, 0);
%
% olddir: directory where old analyses are saved
% newdir: directory where new analyses should be saved
% numreps: number of repetitions - constrained by number of analyses in olddir
% vars2keep: indices of variables which should be kept, i.e. exclude the ones you don't want anymore
% deleteprev: if deleteprev=1 then any results that were already saved in newdir will be deleted. if deleteprev=0 then results in newdir will not be deleted, i.e. analyses will not be rerun

for rep=1:numreps
    for null=1:2
        if null==1
            nullstr='_actual';
        else
            nullstr='_null';
        end
        % load the files from feature thresholding in the original analysis
        load([olddir filesep num2str(rep) nullstr filesep 'design.mat']);
        load([olddir filesep num2str(rep) nullstr filesep 'merit_per_var.mat']);
        % if teh folder already exists from a previous analyses and
        % deletevars was turned on then delete the previous results
        if isequal(deleteprev,1) && exist([newdir filesep num2str(rep) nullstr], 'dir')
            rmdir([newdir filesep num2str(rep) nullstr]);
        end
        % if the folder doesn't already exist make it (it would still exist if deletevars
        % is turned off and part or all of the analysis was already conducted)
        if ~exist([newdir filesep num2str(rep) nullstr], 'dir')
            mkdir([newdir filesep num2str(rep) nullstr]);
        end
        % if deletevars is turned off and part or all of the analysis was
        % already conducted then this bit will not be done
        if ~exist([newdir filesep num2str(rep) nullstr filesep 'merit_per_var.mat'], 'file')
            design.data=design.data(:,vars2keep);
            design.vars=design.vars(vars2keep);
            design.saveto=[newdir filesep num2str(rep) nullstr];
            save([design.saveto filesep 'design.mat'], 'design');
            oldmerit=merit_per_var;
            clear merit_per_var
            for fold=1:(design.numFolds*design.numFolds)
                merit_per_var{fold}=oldmerit{fold}(vars2keep);
            end
            save([design.saveto filesep 'merit_per_var.mat'], 'merit_per_var');
        end
        stats=eliminatereruns([newdir filesep num2str(rep) nullstr]);
    end
end