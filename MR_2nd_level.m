function [LHmerit, vars2use] = MR_2nd_level(design, pass_vars, folds2run)

% latest update: may 30th 2017

cd(design.saveto);

fprintf('Performing Middle Fold Analysis with %d subjects and %d predictors\n', size(design.data));

design.numAlphas=1;
design.numLambdas=1;
tmpLHmerit=cell(design.numFolds*design.numFolds,1);
tmpvars2use=cell(design.numFolds*design.numFolds,1);
for folds=1:(design.numFolds*design.numFolds)
    tmpvars2use{folds}=zeros(design.numFolds, size(design.merit_thresholds,2),design.nvars);
    tmpLHmerit{folds}=NaN(design.numFolds, size(design.merit_thresholds,2));
end
parfor folds=1:length(folds2run)
    [mainfold, subfolds]=ind2sub([design.numFolds design.numFolds], folds2run(folds));
    %identify which threshold constellations are duplicates
    combinations2do=[];
    tofillin=[];
    ROI_duplicate_match=[];
    merit_duplicate_match=[];
    for thresholds=1:design.numFolds*size(design.merit_thresholds,2)
        [ROIcrit meritthresh_looper]=ind2sub([design.numFolds size(design.merit_thresholds,2)], thresholds);
        ROIcrit=ROIcrit-1;  ROIcrit_looper=design.numFolds-ROIcrit;
        if isempty(find(pass_vars(ROIcrit_looper, meritthresh_looper, mainfold, subfolds,:)==1))==0
            [prev_ROIcritlooper, prev_meritthresh]=check_duplicate(design, pass_vars, ROIcrit_looper, meritthresh_looper, mainfold, subfolds);
            if isempty(prev_ROIcritlooper)==1
                %these are the combinations that we need to calculate
                combinations2do=[combinations2do,thresholds];
            else
                %these are combinations that have duplicates so we can
                %fill them in later
                tofillin=[tofillin,thresholds];
                ROI_duplicate_match=[ROI_duplicate_match,prev_ROIcritlooper];
                merit_duplicate_match=[ROI_duplicate_match,prev_meritthresh];
            end
        end
    end
    
    %get training and test set for the subfold
    trainsubjectsTune = find(design.mainfold ~= mainfold & design.subfolds(:,mainfold) ~=subfolds);
    testsubjectsTune = find(design.mainfold ~= mainfold & design.subfolds(:,mainfold)==subfolds);
    
    %loop over all combinations that need to be calculated
    for thresholds=1:length(combinations2do)
        [ROIcrit, meritthresh_looper]=ind2sub([design.numFolds size(design.merit_thresholds,2)], combinations2do(thresholds));
        ROIcrit=ROIcrit-1;  ROIcrit_looper=design.numFolds-ROIcrit;
        %get the feature set
        tmpVars2pick=find(pass_vars(ROIcrit_looper, meritthresh_looper, mainfold, subfolds,:)==1);
        %add the covariates
        Vars2pick=[tmpVars2pick; [size(design.data,2)+1:size(design.data,2)+size(design.extradata,2)]'];
        
        tmpmerit=zeros(design.nboot,1);
        ttmpvars2use=zeros(design.nboot,design.nvars);
        if design.nboot>1
            for bootct=1:design.nboot
                train=0;
                while length(train)<2
                    %get bootstrap data subset
                    [Xboot,Yboot]=bootstrapal([design.data(trainsubjectsTune,Vars2pick), design.extradata(trainsubjectsTune,:)],design.outcome(trainsubjectsTune),design.Ratio);
                    %make sure we have multiple training subjects
                    train=unique(Yboot);
                end
                
                [b,bint,r,rint,stats] = regress(Yboot,Xboot);
                usedvars=zeros(design.numLambdas,design.nvars);
                tmp=zeros(size(design.data,2),1);
                for nchk=1:length(tmpVars2pick)
                    if b(nchk)~=0
                        tmp(Vars2pick(nchk))=1;
                    end
                end
                usedvars=tmp;
                %make outcome predictions
                getprobstune = glmval(squeeze([0;b]),design.data(testsubjectsTune,Vars2pick),design.link);
                
                %get evaluation metric
                tmp1 = -sqrt(abs(design.data(testsubjectsTune)-getprobstune)'*abs(design.data(testsubjectsTune)-getprobstune)/length(design.data(testsubjectsTune)));
                merit=tmp1;
                ttmpvars2use(bootct,:)=usedvars;
                tmpmerit(bootct)=merit;
            end
            %aggregate the results from all bootstrap iterations
            if strcmp(design.bagcrit, 'cdf')==1
                %get cumulative distribution function
                tmpy=cdf('norm', tmpmerit, mean(tmpmerit), std(tmpmerit));
                %find the value at which the likelihood of overestimating the result is equivalent to design.siglevel (normally 0.05)
                findLH=tmpy(find(tmpy<=design.siglevel));
                findLH=max(findLH);
                if isempty(findLH)==0
                    ttmpLHmerit=tmpmerit(find(tmpy==findLH));
                else
                    ttmpLHmerit=NaN;
                end
                tmpLHmerit{folds}(ROIcrit_looper,meritthresh_looper)=ttmpLHmerit(1);
            elseif strcmp(design.bagcrit, 'median')==1
                tmpLHmerit{folds}(ROIcrit_looper,meritthresh_looper)=median(tmpmerit);
            else
                disp('Please enter a valid method for bagging (median or cdf)')
            end
            %find in how many bootstrap iterations the elastic net excluded each variable
            tmpfindvars=sum(squeeze(ttmpvars2use),1);
            if length(tmpfindvars)~=size(design.data,2)
                tmpfindvars=sum(squeeze(ttmpvars2use),2);
            end
            tmpfind=zeros(size(design.data,2),1);
            %get variables that were included in 50% of iterations
            tmpfind(find(tmpfindvars>=design.nboot*.5))=1;
            tmpvars2use{folds}(ROIcrit_looper,meritthresh_looper,:)=logical(tmpfind);
        else
            X2use=[design.data, design.extradata];
            [b,bint,r,rint,stats] = regress(design.outcome(trainsubjectsTune),X2use(trainsubjectsTune, Vars2pick));
            usedvars=zeros(design.numLambdas,design.nvars);
            tmp=zeros(size(design.data,2),1);
            for nchk=1:length(tmpVars2pick)
                if b(nchk)~=0
                    tmp(Vars2pick(nchk))=1;
                end
            end
            usedvars=tmp;
            getprobstune = glmval(squeeze([0;b]),X2use(testsubjectsTune,Vars2pick),design.link);
            tmp1 = -sqrt(abs(design.outcome(testsubjectsTune)-getprobstune)'*abs(design.outcome(testsubjectsTune)-getprobstune)/length(design.outcome(testsubjectsTune)));
            merit=tmp1;
            tmpvars2use{folds}(ROIcrit_looper,meritthresh_looper,:)=logical(usedvars)';
            tmpLHmerit{folds}(ROIcrit_looper,meritthresh_looper)=merit;
        end
    end
disp(['finished mainfold ' num2str(mainfold) ' subfold ' num2str(subfolds)]);
 save_folds_in_parfor(tmpLHmerit{folds}, tmpvars2use{folds}, mainfold, subfolds, design);
end
% vars2use=zeros(design.numFolds, design.numFolds, design.numFolds, size(design.merit_thresholds,2),design.nvars);
% LHmerit=NaN(design.numFolds, design.numFolds, design.numFolds, size(design.merit_thresholds,2));
% for folds=1:(design.numFolds*design.numFolds)
%     [mainfold subfolds]=ind2sub([design.numFolds design.numFolds], folds);
%     LHmerit(mainfold, subfolds,:,:)=tmpLHmerit{folds};
%     vars2use(mainfold,subfolds,:,:,:)=tmpvars2use{folds};
% end
% vars2use=logical(vars2use);
% save([pwd filesep 'LHmerit.mat'], 'LHmerit')
% try
%     save([pwd filesep 'vars2use.mat'], 'vars2use')
% catch
%     save([pwd filesep 'vars2use.mat'], 'vars2use', '-v7.3')
% end
end

function [prev_ROIcritlooper, prev_AUCthresh]=check_duplicate(design, AUCs2pick, ROIcrit_looper, AUCthresh_looper, mainfold, subfolds)
prev_ROIcritlooper=[];
prev_AUCthresh=[];

a=ROIcrit_looper;
b=AUCthresh_looper;
c=mainfold;
d=subfolds;

if isempty(find(AUCs2pick(a,b,c,d,:)==1))==0
    for e=(a+1):design.numFolds %search one fold up
        if isempty(prev_ROIcritlooper)==1
            for f=1:length(design.merit_thresholds)%AUCthreshold
                if isequal(find(AUCs2pick(a,b,c,d,:)==1), find(AUCs2pick(e,f,c,d,:)==1))==1
                    prev_ROIcritlooper=e;
                    prev_AUCthresh=f;
                end
            end
        end
    end
end
end


function save_folds_in_parfor(tmpLHmerit, tmpvars2use,mainfold, subfolds, design)
%rather than holding all data in memory files are saved intermittently.
% see eliminatereruns.m for code to use these intermittently saved files to
% restart an analysis if it has been halted/killed/crashed
cd(design.saveto);
    LHmerit2save=tmpLHmerit;
    vars2use2save=tmpvars2use;
    save([design.saveto filesep 'Second_level_results_' num2str(mainfold) '_' num2str(subfolds)], 'LHmerit2save', 'vars2use2save');
end
