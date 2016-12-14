function [design, pass_vars]=RAFT_do_thresh(design, merit_per_var)
% thresholding step of Regularized Adaptive Feature Thresholding
%
% for comments and questions please contact lee.jollans@gmail.com

fprintf('Determining thresholds with %d subjects and %d predictors\n', size(design.data));

pass_vars=zeros(design.numFolds, 10, design.numFolds, design.numFolds, design.nvars);

for vars=1:design.nvars
    for folds=1:design.numFolds*design.numFolds
        [mainfold subfolds]=ind2sub([design.numFolds design.numFolds], folds);
        featuremerit(mainfold, subfolds, vars)=merit_per_var{folds}(vars);
    end
end

for mainfold=1:design.numFolds
    maxcollect=[];
% strictest threshold -- at least one feature left in each fold
for subfolds=1:design.numFolds
    maxcollect=[maxcollect, max(featuremerit(mainfold, subfolds, :))];
end
maxthresh=min(maxcollect)-.001;
% most liberal threshold -- at least one feature that appears in all subfolds in a mainfold, for all mainfolds
    tmp=squeeze(featuremerit(mainfold,:,:));
    for n=1:design.nvars
    featuremax(n)=min(tmp(:,n));
    end
    minthresh=max(featuremax)-.001;
if minthresh==maxthresh
    minthresh=minthresh-.01;
end
    design.merit_thresholds(mainfold,:)=linspace(minthresh, maxthresh, 10); 
    for threshold=1:10 %number of 'merit' thresholds hardcoded as 10
        tmppass_vars=zeros(design.numFolds, design.nvars);
        for subfolds=1:design.numFolds
            track=squeeze(featuremerit(mainfold,subfolds,:))-design.merit_thresholds(mainfold,threshold);
            tmppass_vars(subfolds, find(track>0))=1;
        end
            occurrence_chk=squeeze(sum(tmppass_vars,1));
            for occurrence=1:design.numFolds
                idx=find(occurrence_chk>=occurrence);
                for subfold=1:design.numFolds
                    pass_vars(occurrence, threshold, mainfold, subfold, idx)=tmppass_vars(subfold,idx);
                end
            end
        end
    end
    cd(design.saveto)
    save('pass_vars', 'pass_vars')
    save('design', 'design')
end
