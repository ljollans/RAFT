function RAFT_2nd_Level(design, pass_vars, folds2run)
% Second level for Regularized Adaptive Feature Thresholding
% Here models for each combination of parameters and thresholds are
% evaluated and the superior models are selected at the third level
%
% for comments and questions please contact lee.jollans@gmail.com

cd(design.saveto);
design.nvars=size(design.data,2);
design.extradata(find(isnan(design.extradata)==1))=0;

fprintf('Performing Model Evaluations with %d subjects and %d predictors\n', size(design.data));

parfor folds=1:(design.numFolds*design.numFolds)
    if isempty(find(folds==folds2run))==0;
        
            tmpLHmerit=NaN(design.numFolds, size(design.merit_thresholds,2),design.numAlphas, design.numLambdas);
    tmpvars2use=zeros(design.numFolds, size(design.merit_thresholds,2),design.numAlphas, design.numLambdas,size(design.data,2));
    tmplambda_values=NaN(design.numFolds, size(design.merit_thresholds,2),design.numAlphas, design.numLambdas);
        
    [mainfold, subfolds]=ind2sub([design.numFolds design.numFolds], folds);
    %identify which threshold constellations are duplicates
    combinations2do=[];
    tofillin=[];
    ROI_duplicate_match=[];
    merit_duplicate_match=[];
    for thresholds=1:design.numFolds*size(design.merit_thresholds,2)
        [ROIcrit meritthresh_looper]=ind2sub([design.numFolds size(design.merit_thresholds,2)], thresholds);
        if isempty(find(pass_vars(ROIcrit, meritthresh_looper, mainfold, subfolds,:)==1))==0
            [prev_ROIcritlooper, prev_meritthresh]=check_duplicate(design, pass_vars, ROIcrit, meritthresh_looper, mainfold, subfolds);
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
        %get the feature set
        tmpVars2pick=find(pass_vars(ROIcrit, meritthresh_looper, mainfold, subfolds,:)==1);
        %add the covariates
        Vars2pick=[tmpVars2pick; [size(design.data,2)+1:size(design.data,2)+size(design.extradata,2)]'];
        
        %loop over alpha values
        for alpha_looper=1:length(design.alpha)
            ttmpvars2use=zeros(design.nboot,design.numLambdas,size(design.data,2));
            tmpmerit=NaN(design.nboot,design.numLambdas);
            tmplambda2use=zeros(design.nboot,design.numLambdas);
            %bootstrap iterations
            if design.nboot>1
                for bootct=1:design.nboot
                    %glmnet options
                    options=glmnetSet;
                    options.standardize=true;
                    options.nlambda=length(design.lambda);
                    options.lambda=design.lambda;
                    options.alpha=design.alpha(alpha_looper);
                    train=0;
                    while train==0
                        %get bootstrap data subset
                        try
                        [Xboot,Yboot]=bootstrapal([design.data([trainsubjectsTune; testsubjectsTune],:), design.extradata([trainsubjectsTune; testsubjectsTune],:)],design.outcome([trainsubjectsTune; testsubjectsTune]),design.Ratio);
                        catch ME
                            [Xboot,Yboot]=bootstrapal([design.data([trainsubjectsTune; testsubjectsTune],:), design.extradata([trainsubjectsTune; testsubjectsTune],:)],design.outcome([trainsubjectsTune; testsubjectsTune])',design.Ratio);
                        end
                        %make sure we have multiple training subjects
                        if length(unique(Yboot(length(trainsubjectsTune)+1:end)))>1 
                            chk1=1;
                        else
                            chk1=0;
                        end
                        %make sure we have multiple levels in at least one
                        %predictor
                        t1=[];
                        for v=1:length(Vars2pick)
                            t1(v)=length(unique(Xboot(1:length(trainsubjectsTune), Vars2pick(v))));
                        end
                        if isempty(find(t1>1))==0
                            chk2=1;
                        else
                            chk2=0;
                        end
                        if chk1==1 && chk2==1
                            train=1;
                        end
                    end
                    
                    %fit the elastic net and get beta and intercept values
                    try
                    fit=glmnet(Xboot(1:length(trainsubjectsTune), Vars2pick),Yboot(1:length(trainsubjectsTune)),design.family,options);
                    catch
                        pause
                    end
                    B0=fit.beta; intercept=fit.a0;
                    
                    switch(design.type)
                        case 'linear',
                            %if initial lambda value was too high, adjust it
                            while size(B0,2)<options.nlambda
                                options.lambda=linspace(min(options.lambda),max(fit.lambda)+(max(fit.lambda)/5),design.numLambdas);
                                fit=glmnet(Xboot(1:length(trainsubjectsTune), Vars2pick),Yboot(1:length(trainsubjectsTune)),design.family,options);
                                B0=fit.beta; intercept=fit.a0;
                            end
                        case 'logistic',
                            %if initial lambda value was too high, adjust it
                            while size(B0,2)<options.nlambda
                                options.lambda=linspace(min(options.lambda),max(fit.lambda)+(max(fit.lambda)/5),design.numLambdas);
                                fit=glmnet(Xboot(1:length(trainsubjectsTune), Vars2pick),Yboot(1:length(trainsubjectsTune)),design.family,options);
                                B0=fit.beta; intercept=fit.a0;
                            end
                    end

                    tmplambda2use(bootct,:)=fit.lambda;
                    try b=[intercept'; B0];
                    catch
                        b=[intercept; B0];
                    end
                    merit=NaN(1,design.numLambdas);
                    usedvars=zeros(design.numLambdas,size(design.data,2));
                    
                    for lambda_looper=1:length(design.lambda)
                        bb=b(:,lambda_looper);
                        %identify what variables the elastic net excluded
                        tmp=zeros(size(design.data,2),1);
                        for nchk=1:length(tmpVars2pick)
                            if bb(nchk+1)~=0 %+1 because of the intercept
                                tmp(Vars2pick(nchk))=1;
                            end
                        end
                        usedvars(lambda_looper,:)=tmp;
                        if sum(tmp)==0
                            tmp(Vars2pick)=1;
                        end
                        
                        %make outcome predictions
                        getprobstune = glmval(squeeze(bb),Xboot(length(trainsubjectsTune)+1:end,Vars2pick),design.link);
                        
                        %get evaluation metric
                        switch(design.type)
                            case 'linear',
                                truth=Yboot(length(trainsubjectsTune)+1:end);
                                tmp1 = -sqrt(abs(truth-getprobstune)'*abs(truth-getprobstune)/length(truth));
                                merit(lambda_looper)=tmp1;
                            case 'logistic',
                                switch(design.balanced)
                                    case 'balanced'
                                        [tmp1] = fastAUC(logical(Yboot(length(trainsubjectsTune)+1:end)),getprobstune,0);
                                        merit(lambda_looper)=tmp1;
                                    case 'unbalanced'
                                        [prec, tpr] = prec_rec_rob_mod(getprobstune, Yboot(length(trainsubjectsTune)+1:end),'tstPrecRec', 'plotPR',0, 'numThresh',100);
                                        fscore=(prec.*tpr)./(prec+tpr);
                                        merit(lambda_looper)=max(fscore);
                                end
                        end
                    end
                    ttmpvars2use(bootct,:,:)=usedvars;
                    tmpmerit(bootct,:)=merit;
                end
                %collect the lambda values which the elastic net used
                tmplambda_values(ROIcrit, meritthresh_looper, alpha_looper,:)=mean(tmplambda2use);
                %aggregate the results from all bootstrap iterations
                for lambda_looper=1:length(design.lambda)
                    if strcmp(design.bagcrit, 'cdf')==1
                        %get cumulative distribution function
                        tmpy=cdf('norm', tmpmerit(:,lambda_looper), mean(tmpmerit(:,lambda_looper)), std(tmpmerit(:,lambda_looper)));
                        %find the value at which the likelihood of overestimating the result is equivalent to design.siglevel (normally 0.05)
                        findLH=tmpy(find(tmpy<=design.siglevel));
                        findLH=max(findLH);
                        if isempty(findLH)==0
                            ttmpLHmerit=tmpmerit(find(tmpy==findLH),lambda_looper);
                        else
                            ttmpLHmerit=NaN;
                        end
                        tmpLHmerit(ROIcrit, meritthresh_looper, alpha_looper,lambda_looper)=ttmpLHmerit(1);
                    elseif strcmp(design.bagcrit, 'median')==1
                        tmpLHmerit(ROIcrit, meritthresh_looper, alpha_looper,lambda_looper)=median(tmpmerit(:,lambda_looper));%CHANGED TO MEDIAN 25TH APRIL!!
                    else
                        disp('Please enter a valid method for bagging (median or cdf)')
                    end
                    
                    %find in how many bootstrap iterations the elastic net excluded each variable
                    tmpfindvars=sum(squeeze(ttmpvars2use(:,lambda_looper,:)),1);
                    if length(tmpfindvars)~=size(design.data,2)
                        tmpfindvars=sum(squeeze(ttmpvars2use(:,lambda_looper,:)),2);
                    end
                    tmpfind=zeros(size(design.data,2),1);
                    %get variables that were included in more than 50% of bags
                    tmpfind(find(tmpfindvars>design.nboot*.5))=1;
                    tmpvars2use(ROIcrit, meritthresh_looper, alpha_looper,lambda_looper,:)=logical(tmpfind);
                end
            else
                %glmnet options
                options=glmnetSet;
                options.standardize=true;
                options.nlambda=length(design.lambda);
                options.lambda=design.lambda;
                options.alpha=design.alpha(alpha_looper);
                X2use=[design.data, design.extradata];
                %fit the elastic net and get beta and intercept values
                fit=glmnet(X2use(trainsubjectsTune, Vars2pick),design.outcome(trainsubjectsTune),design.family,options);
                B0=fit.beta; intercept=fit.a0;
                
                switch(design.type)
                    case 'linear',
                        %if initial lambda value was too high, adjust it
                        while size(B0,2)<options.nlambda
                            options.lambda=linspace(min(options.lambda),max(fit.lambda)+(max(fit.lambda)/5),design.numLambdas);
                            fit=glmnet(X2use(trainsubjectsTune, Vars2pick),design.outcome(trainsubjectsTune),design.family,options);
                            B0=fit.beta; intercept=fit.a0;
                        end
                    case 'logistic',
                        %if initial lambda value was too high, adjust it
                        while size(B0,2)<options.nlambda
                            options.lambda=linspace(min(options.lambda),max(fit.lambda)+(max(fit.lambda)/5),design.numLambdas);
                            fit=glmnet(X2use(trainsubjectsTune, Vars2pick),design.outcome(trainsubjectsTune),design.family,options);
                            B0=fit.beta; intercept=fit.a0;
                        end
                end
                tmplambda_values(ROIcrit, meritthresh_looper, alpha_looper,:)=fit.lambda;
                try b=[intercept'; B0];
                catch
                    b=[intercept; B0];
                end
                
                for lambda_looper=1:length(design.lambda)
                    merit=zeros(1,design.numLambdas);
                    usedvars=zeros(design.numLambdas,size(design.data,2));
                    bb=b(:,lambda_looper);
                    
                    %identify what variables the elastic net excluded
                    tmp=zeros(size(design.data,2),1);
                    for nchk=1:length(tmpVars2pick)
                        if bb(nchk+1)~=0
                            tmp(Vars2pick(nchk))=1;
                        end
                    end
                    tmpvars2use(ROIcrit, meritthresh_looper, alpha_looper,lambda_looper,:)=logical(tmp);
                    
                    %make outcome predictions
                    getprobstune = glmval(squeeze(bb),X2use(testsubjectsTune,Vars2pick),design.link);
                    
                    %get evaluation metric
                    switch(design.type)
                        case 'linear',
                            tmp1 = -sqrt(abs(design.outcome(testsubjectsTune)-getprobstune)'*abs(design.outcome(testsubjectsTune)-getprobstune)/length(design.outcome(testsubjectsTune)));
                            merit(lambda_looper)=tmp1;
                        case 'logistic',
                            switch(design.balanced)
                                case 'balanced'
                                    [tmp1] = fastAUC(logical(design.outcome(testsubjectsTune)),getprobstune,0);
                                    merit(lambda_looper)=tmp1;
                                case 'unbalanced'
                                    [prec, tpr] = prec_rec_rob_mod(getprobstune, design.outcome(testsubjectsTune),'tstPrecRec', 'plotPR',0, 'numThresh',100);
                                    fscore=(prec.*tpr)./(prec+tpr);
                                    merit(lambda_looper)=max(fscore);
                            end
                    end
                end
                %                 end
                tmpLHmerit(ROIcrit, meritthresh_looper, alpha_looper,:)=merit;
            end
        end
        %collect the lambda values which the elastic net used
        tmplambda_values(ROIcrit, meritthresh_looper, alpha_looper,:)=mean(tmplambda2use);
    end
    
    %% fill in duplicates
    for filler=1:length(tofillin)
        [ROIcrit meritthresh_looper]=ind2sub([design.numFolds size(design.merit_thresholds,2)], tofillin(filler));
        tmplambda_values(ROIcrit, meritthresh_looper, :,:)=tmplambda_values(ROI_duplicate_match(filler), merit_duplicate_match(filler), :,:);
        tmpLHmerit(ROIcrit, meritthresh_looper, :,:)=tmpLHmerit(ROI_duplicate_match(filler), merit_duplicate_match(filler), :,:);
        tmpvars2use(ROIcrit, meritthresh_looper, :,:,:)=tmpvars2use(ROI_duplicate_match(filler), merit_duplicate_match(filler), :,:,:);
    end
    
    disp(['finished mainfold ' num2str(mainfold) ' subfold ' num2str(subfolds)]);
    save_folds_in_parfor(tmpLHmerit, tmpvars2use, tmplambda_values, mainfold, subfolds, design);
    end
end
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

function save_folds_in_parfor(tmpLHmerit, tmpvars2use, tmplambda_values, mainfold, subfolds, design)
%rather than holding all data in memory files are saved intermittently.
% see eliminatereruns.m for code to use these intermittently saved files to
% restart an analysis if it has been halted/killed/crashed
cd(design.saveto);
    LHmerit2save=tmpLHmerit;
    vars2use2save=tmpvars2use;
    lambda_values2save=tmplambda_values;
    try
    save([design.saveto filesep 'Second_level_results_' num2str(mainfold) '_' num2str(subfolds)], 'LHmerit2save', 'vars2use2save', 'lambda_values2save');
    catch ME
        try
            mkdir([design.saveto filesep 'tmp']);
            save([design.saveto filesep 'tmp' filesep 'Second_level_results_' num2str(mainfold) '_' num2str(subfolds)], 'LHmerit2save', 'vars2use2save', 'lambda_values2save');
        catch ME
            disp('I do not know what this saving error is! help!')
        end
    end
end

function [prec, tpr, fpr, thresh] = prec_rec_rob_mod(score, target,titleofplot, varargin)
% PREC_REC - Compute and plot precision/recall and ROC curves.
%
%   PREC_REC(SCORE,TARGET), where SCORE and TARGET are equal-sized vectors,
%   and TARGET is binary, plots the corresponding precision-recall graph
%   and the ROC curve.
%
%   Several options of the form PREC_REC(...,'OPTION_NAME', OPTION_VALUE)
%   can be used to modify the default behavior.
%      - 'instanceCount': Usually it is assumed that one line in the input
%                         data corresponds to a single sample. However, it
%                         might be the case that there are a total of N
%                         instances with the same SCORE, out of which
%                         TARGET are classified as positive, and (N -
%                         TARGET) are classified as negative. Instead of
%                         using repeated samples with the same SCORE, we
%                         can summarize these observations by means of this
%                         option. Thus it requires a vector of the same
%                         size as TARGET.
%      - 'numThresh'    : Specify the (maximum) number of score intervals.
%                         Generally, splits are made such that each
%                         interval contains about the same number of sample
%                         lines.
%      - 'holdFigure'   : [0,1] draw into the current figure, instead of
%                         creating a new one.
%      - 'style'        : Style specification for plot command.
%      - 'plotROC'      : [0,1] Explicitly specify if ROC curve should be
%                         plotted.
%      - 'plotPR'       : [0,1] Explicitly specify if precision-recall curve
%                         should be plotted.
%      - 'plotBaseline' : [0,1] Plot a baseline of the random classifier.
%
%   By default, when output arguments are specified, as in
%         [PREC, TPR, FPR, THRESH] = PREC_REC(...),
%   no plot is generated. The arguments are the score thresholds, along
%   with the respective precisions, true-positive, and false-positive
%   rates.
%
%   Example:
%
% x1 = rand(1000, 1);
% y1 = round(x1 + 0.5*(rand(1000,1) - 0.5));
% prec_rec(x1, y1);
% x2 = rand(1000,1);
% y2 = round(x2 + 0.75 * (rand(1000,1)-0.5));
% prec_rec(x2, y2, 'holdFigure', 1);
% legend('baseline','x1/y1','x2/y2','Location','SouthEast');

% Copyright ï¿½ 9/22/2010 Stefan Schroedl
% Updated     3/16/2010

optargin = size(varargin, 2);
stdargin = nargin - optargin;

if stdargin < 3
    error('at least 2 arguments required');
end

% parse optional arguments
num_thresh = -1;
hold_fig = 0;
plot_roc = (nargout <= 0);
plot_pr  = (nargout <= 0);
instance_count = -1;
style = '';
plot_baseline = 1;

i = 1;
while (i <= optargin)
    if (strcmp(varargin{i}, 'numThresh'))
        if (i >= optargin)
            error('argument required for %s', varargin{i});
        else
            num_thresh = varargin{i+1};
            i = i + 2;
        end
    elseif (strcmp(varargin{i}, 'style'))
        if (i >= optargin)
            error('argument required for %s', varargin{i});
        else
            style = varargin{i+1};
            i = i + 2;
        end
    elseif (strcmp(varargin{i}, 'instanceCount'))
        if (i >= optargin)
            error('argument required for %s', varargin{i});
        else
            instance_count = varargin{i+1};
            i = i + 2;
        end
    elseif (strcmp(varargin{i}, 'holdFigure'))
        if (i >= optargin)
            error('argument required for %s', varargin{i});
        else
            if ~isempty(get(0,'CurrentFigure'))
                hold_fig = varargin{i+1};
            end
            i = i + 2;
        end
    elseif (strcmp(varargin{i}, 'plotROC'))
        if (i >= optargin)
            error('argument required for %s', varargin{i});
        else
            plot_roc = varargin{i+1};
            i = i + 2;
        end
    elseif (strcmp(varargin{i}, 'plotPR'))
        if (i >= optargin)
            error('argument required for %s', varargin{i});
        else
            plot_pr = varargin{i+1};
            i = i + 2;
        end
    elseif (strcmp(varargin{i}, 'plotBaseline'))
        if (i >= optargin)
            error('argument required for %s', varargin{i});
        else
            plot_baseline = varargin{i+1};
            i = i + 2;
        end
    elseif (~ischar(varargin{i}))
        error('only two numeric arguments required');
    else
        error('unknown option: %s', varargin{i});
    end
end

[nx,ny]=size(score);

if (nx~=1 && ny~=1)
    error('first argument must be a vector');
end

[mx,my]=size(target);
if (mx~=1 && my~=1)
    error('second argument must be a vector');
end

score  =  score(:);
target = target(:);

if (length(target) ~= length(score))
    error('score and target must have same length');
end

if (instance_count == -1)
    % set default for total instances
    instance_count = ones(length(score),1);
    target = max(min(target(:),1),0); % ensure binary target
else
    if numel(instance_count)==1
        % scalar
        instance_count = instance_count * ones(length(target), 1);
    end
    [px,py] = size(instance_count);
    if (px~=1 && py~=1)
        error('instance count must be a vector');
    end
    instance_count = instance_count(:);
    if (length(target) ~= length(instance_count))
        error('instance count must have same length as target');
    end
    target = min(instance_count, target);
end

if num_thresh < 0
    % set default for number of thresholds
    score_uniq = unique(score);
    num_thresh = min(length(score_uniq), 100);
    if num_thresh<=1
        num_thresh=3;
    end
end

qvals = (1:(num_thresh-1))/num_thresh;
thresh = [min(score) quantile(score,qvals)];
% remove identical bins
thresh = sort(unique(thresh),2,'descend');
total_target = sum(target);
total_neg = sum(instance_count - target);

prec = zeros(length(thresh),1);
tpr  = zeros(length(thresh),1);
fpr  = zeros(length(thresh),1);
for i = 1:length(thresh)
    idx     = (score >= thresh(i));
    fpr(i)  = sum(instance_count(idx) - target(idx));
    tpr(i)  = sum(target(idx)) / total_target;
    prec(i) = sum(target(idx)) / sum(instance_count(idx));
    
end
fpr = fpr / total_neg;

end
