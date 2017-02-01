function [params2pick, Merit, Beta, Vars2pick_main,  GetProbs, design, stats] = RAFT_Model_Selection(design, vars2use, LHmerit, lambda_values, pass_vars, use_pass_vars, startROIthresh, startmeritthresh)
% select the best performing model across subfolds for each main CV fold
% and applies it to the test set
%
% for comments and questions please contact lee.jollans@gmail.com

design.prediction=[];
fprintf('Performing Mainfold Analysis with %d subjects and %d predictors\n', size(design.data));
design.nvars=size(design.data,2);
outcome2use=design.outcome;

options=glmnetSet;
options.standardize=true;
options.nlambda=1;

[params2pick paramidx subfoldidx lambdavalues] = findbestparams(LHmerit, design, lambda_values, startROIthresh, startmeritthresh);
intalpha=4; intlambda=3; intmerit=2; introi=1;
params2use=params2pick;
for mainfold=1:design.numFolds
    %%
    Vars2pick_main{mainfold}=[];
    chkr=0;
    
    paramidx_perfold=squeeze(subfoldidx(:,mainfold,:));
    [rf, cf]=find(isnan(paramidx_perfold)==1);
    for n=1:length(cf)
        paramidx_perfold(rf(n), cf(n))=round(nanmean(squeeze(paramidx_perfold(rf(n),:))));
    end
    for subfold=1:design.numFolds
        tmpidx=squeeze(paramidx_perfold(:,subfold));
        folds=sub2ind([design.numFolds design.numFolds], mainfold, subfold);
        try
            if use_pass_vars==0
                trackvars(mainfold,subfold,:)=squeeze(vars2use(mainfold, subfold,tmpidx(introi),tmpidx(intmerit),tmpidx(intalpha),tmpidx(intlambda),:));
            elseif use_pass_vars==1
                trackvars(mainfold,subfold,:)=squeeze(pass_vars(tmpidx(introi),tmpidx(intmerit),mainfold, subfold,:));
            end
        catch ME
            disp(mainfold)
        end
    end
    varsscount(mainfold,:)=sum(trackvars(mainfold,:,:),2);
    Vars2pick_main{mainfold}=find(varsscount(mainfold,:)>=params2use(introi,mainfold));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % in case no variables were found matching the identified
    % parameters
    if length(Vars2pick_main{mainfold})<1
        
        %% make an array identifying for what parameter combinations
        %there would be variables to use
        pickr=zeros(design.numFolds, size(design.merit_thresholds,2), design.numAlphas, design.numLambdas);
        while sum(sum(sum(sum(pickr))))==0
            for idxs=1:[design.numFolds*size(design.merit_thresholds,2)*design.numAlphas*design.numLambdas*design.numFolds]
                [roiidxs meritidxs alphaidxs lambdaidxs subfold]=ind2sub([design.numFolds size(design.merit_thresholds,2) design.numAlphas design.numLambdas design.numFolds], idxs);
                folds=sub2ind([design.numFolds design.numFolds], mainfold, subfold);
                if use_pass_vars~=1
                    tmptrackvars(roiidxs, meritidxs, alphaidxs, lambdaidxs,subfold,:)=squeeze(vars2use(mainfold,subfold,roiidxs, meritidxs, alphaidxs, lambdaidxs,:));
                else
                    tmptrackvars(roiidxs, meritidxs, alphaidxs, lambdaidxs,subfold,:)=squeeze(pass_vars(roiidxs, meritidxs, mainfold, subfold,:));
                end
            end
            %take out the parameter combination that was originally
            %tried
            for subfold=1:design.numFolds
                tmpidx=squeeze(subfoldidx(:,mainfold, subfold));
                tmptrackvars(tmpidx(introi),tmpidx(intmerit),tmpidx(intalpha),tmpidx(intlambda),subfold,:)=zeros(size(design.vars));
            end
            for idxs=1:[design.numFolds*size(design.merit_thresholds,2)*design.numAlphas*design.numLambdas*design.numFolds]
                [roiidxs meritidxs alphaidxs lambdaidxs subfold]=ind2sub([design.numFolds size(design.merit_thresholds,2) design.numAlphas design.numLambdas design.numFolds], idxs);
                varscount=sum(tmptrackvars(roiidxs, meritidxs, alphaidxs, lambdaidxs,subfold,:));
                pickr(roiidxs, meritidxs, alphaidxs, lambdaidxs)=not(isempty(find(varscount>=roiidxs)));
            end
            if sum(sum(sum(sum(pickr))))==0 && use_pass_vars==0
                use_pass_vars=1;
            elseif sum(sum(sum(sum(pickr))))==0 && use_pass_vars==1
                error('apparently there was only one parameter combination that would have resulted in a non-empty Vars2pick matrix for this mainfold. This should not be possible!')
            end
        end
        while length(Vars2pick_main{mainfold})<1
            %get all the param combinations that have variables
            [roiidx meritidx alphaidx lambdaidx]=ind2sub([design.numFolds size(design.merit_thresholds,2) design.numAlphas design.numLambdas],find(pickr));
            %note in a matrix where these params correspond to the best
            %params for EACH SUBFOLD (not for the mainfold overall!)
            tmpoptimalchk=zeros(design.numFolds,length(roiidx),4);
            for subfold=1:design.numFolds
                tmpoptimalchk(subfold,find(roiidx==paramidx(introi,mainfold,subfold)), introi)=1;
                tmpoptimalchk(subfold,find(meritidx==paramidx(intmerit,mainfold,subfold)), intmerit)=1;
                tmpoptimalchk(subfold,find(lambdaidx==paramidx(intlambda,mainfold,subfold)), intlambda)=1;
                tmpoptimalchk(subfold,find(alphaidx==paramidx(intalpha,mainfold,subfold)), intalpha)=1;
            end
            optimalchk=squeeze(sum(sum(tmpoptimalchk,3)));
            %find maximum overlap
            tt=find(optimalchk==max(optimalchk));
            %for combinations with maximum overlap check which deviates
            %least from MAINFOLD ideal params
            clear diff_params tmpdiff_params
            for diffloop=1:length(tt)
                for subfold=1:design.numFolds
                    tmpdiff_params(diffloop,subfold)=sum(abs(...
                        [roiidx(tt(diffloop)) meritidx(tt(diffloop)) lambdaidx(tt(diffloop)) alphaidx(tt(diffloop))]...
                        -subfoldidx(:,mainfold,subfold)'));
                end
                diff_params(diffloop)=sum(tmpdiff_params(diffloop,:));
            end
            paramidxs=tt(find(diff_params==min(diff_params)));
            %plug in the bestfitting parameter combination
            paramidx2use(introi,mainfold)=roiidx(paramidxs(1));
            paramidx2use(intmerit,mainfold)=meritidx(paramidxs(1));
            paramidx2use(intalpha,mainfold)=alphaidx(paramidxs(1));
            paramidx2use(intlambda,mainfold)=lambdaidx(paramidxs(1));
            for subfold=1:design.numFolds
                lambdacollect(subfold)=lambdavalues(roiidx(paramidxs(1)),meritidx(paramidxs(1)),alphaidx(paramidxs(1)),lambdaidx(paramidxs(1)),mainfold,subfold);
            end
            params2use(introi,mainfold)=roiidx(paramidxs(1));
            params2use(intmerit,mainfold)=design.merit_thresholds(mainfold,meritidx(paramidxs(1)));
            params2use(intalpha,mainfold)=design.alpha(alphaidx(paramidxs(1)));
            params2use(intlambda,mainfold)=mean(lambdacollect);
            
            %check if we have anything in Vars2pick now
            for subfold=1:design.numFolds
                tmpidx=squeeze(paramidx2use(:,mainfold));
                folds=sub2ind([design.numFolds design.numFolds], mainfold, subfold);
                try
                    if use_pass_vars==0
                        trackvars(mainfold,subfold,:)=squeeze(vars2use(mainfold, subfold,tmpidx(introi),tmpidx(intmerit),tmpidx(intalpha),tmpidx(intlambda),:));
                    elseif use_pass_vars==1
                        trackvars(mainfold,subfold,:)=squeeze(pass_vars(tmpidx(introi),tmpidx(intmerit),mainfold, subfold,:));
                    end
                catch ME
                    disp(mainfold)
                end
            end
            varsscount(mainfold,:)=sum(trackvars(mainfold,:,:),2);
            Vars2pick_main{mainfold}=find(varsscount(mainfold,:)>=params2use(introi,mainfold));
            %if for whatever reason we didn't find any variable for
            %Vars2pick take that parameter combination out and try again
            if length(Vars2pick_main{mainfold})<1
                pickr(tmpidx(introi),tmpidx(intmerit),tmpidx(intalpha),tmpidx(intlambda))=0;
            end
        end
    end
    
    %%
    trainsubjects = find(design.mainfold ~= mainfold);
    testsubjects=find(design.mainfold == mainfold);
    
    options.alpha=params2use(intalpha,mainfold);
    options.lambda=params2use(intlambda,mainfold);
    try
        fit=glmnet([design.data(trainsubjects,Vars2pick_main{mainfold}), design.extradata(trainsubjects,:)],design.outcome(trainsubjects),design.family,options);
    catch
        design.outcome=design.outcome';
        fit=glmnet([design.data(trainsubjects,Vars2pick_main{mainfold}), design.extradata(trainsubjects,:)],design.outcome(trainsubjects),design.family,options);
    end
    B0=fit.beta;
    
    switch(design.type)
        case 'linear',
            Beta{mainfold}=[fit.a0; B0];
            GetProbs{mainfold} = glmval(Beta{mainfold},[design.data(testsubjects,Vars2pick_main{mainfold}), design.extradata(testsubjects,:)],design.link);
            design.prediction(testsubjects)=GetProbs{mainfold};
            Merit(mainfold)= -sqrt(abs(design.outcome(testsubjects)-GetProbs{mainfold})'*abs(design.outcome(testsubjects)-GetProbs{mainfold})/length(design.outcome(testsubjects)));
            [r(mainfold), p(mainfold)] = corr(GetProbs{mainfold}, design.outcome(testsubjects));
        case 'logistic',
            Beta{mainfold}=[fit.a0; B0];
            GetProbs{mainfold} = glmval(Beta{mainfold},[design.data(testsubjects,Vars2pick_main{mainfold}), design.extradata(testsubjects,:)],design.link);
            [r p]=corr(design.outcome(testsubjects), GetProbs{mainfold});
            design.prediction(testsubjects)=GetProbs{mainfold};
            switch(design.balanced)
                case 'balanced'
                    [Merit.aucbalanced(mainfold),fpr,tpr] = fastAUC(logical(design.outcome(testsubjects)),GetProbs{mainfold},0);
                case 'unbalanced'
                    try
                        [prec, tpr, fpr, thresh] = prec_rec_rob_mod(GetProbs{mainfold}, logical(design.outcome(testsubjects)),'tstPrecRec', 'plotPR',0);
                        [Merit.AUC_unbalanced(mainfold),fpr2,tpr2] = fastAUC(logical(design.outcome(testsubjects)),GetProbs{mainfold},0);
                    catch ME
                        length(testsubjects);
                    end
                    fscore=(prec.*tpr)./(prec+tpr);
                    Merit.F1score(mainfold)=max(fscore);
            end
    end
end
    cd(design.saveto);
    clear stats
    switch(design.type)
        case 'linear',
            [stats.r, stats.p]=corr(design.prediction', design.outcome);
            stats.mse= -sqrt(abs(design.outcome-design.prediction')'*abs(design.outcome-design.prediction')/length(design.outcome));
            save('Results',  'Merit', 'GetProbs', 'Beta', 'params2pick', 'params2use', 'stats', 'Vars2pick_main', 'design');
        case 'logistic'
            [stats.overall_AUC,stats.fpr,stats.tpr] = fastAUC(logical(design.outcome),design.prediction',0);
            [stats.prec, stats.tpr, stats.fpr, stats.thresh] = prec_rec_rob_mod(design.prediction', logical(design.outcome),'tstPrecRec', 'plotPR',0);
            fscore=(stats.prec.*stats.tpr)./(stats.prec+stats.tpr);
            stats.F1score=max(fscore);
            save('Results',  'Merit', 'GetProbs', 'Beta', 'params2pick', 'params2use', 'stats', 'Vars2pick_main', 'design');
    end
end

function [optim_mainfold paramidx_persubfold subfoldidx lambdamatrix] = findbestparams(AUCTune_in, design, lambda_values, startROIthresh, startAUCthresh)

for folds=1:design.numFolds*design.numFolds
    [mainfold subfolds]=ind2sub([design.numFolds design.numFolds], folds);
    AUCTune(:,:,:,:,mainfold,subfolds)=squeeze(AUCTune_in(mainfold, subfolds, startROIthresh:end, startAUCthresh:end, :,:));
    lambdamatrix(:,:,:,:,mainfold,subfolds)=squeeze(lambda_values(mainfold, subfolds, startROIthresh:end, startAUCthresh:end, :,:));
end

alphas=design.alpha;
ROIthresholds=startROIthresh:design.numFolds;
for mainfold=1:design.numFolds
    AUCthresholds=design.merit_thresholds(mainfold,startAUCthresh:end);
    for subfold=1:design.numFolds
        for roi=1:size(AUCTune,1)
            for auc=1:size(AUCTune,2)
                for alph=1:size(AUCTune,3)
                    for lam=1:size(AUCTune,4)
                        alphamatrix(roi, auc, alph, lam, mainfold, subfold)=alphas(alph);
                        AUCthreshmatrix(roi,auc,alph,lam, mainfold, subfold)=AUCthresholds(auc);
                        ROIcritmatrix(roi,auc,alph,lam, mainfold, subfold)=ROIthresholds(roi);
                    end
                end
            end
        end
    end
end
win=zeros(4,design.numFolds);

%% find the best parameters per subfold

for mainfold=1:design.numFolds
    %%
    winsub=[];
    for subfold=1:design.numFolds
        critvarperfold=squeeze(AUCTune(:,:,:,:,mainfold,subfold));
        lambdaperfold=squeeze(lambdamatrix(:,:,:,:,mainfold,subfold));
        alphaperfold=squeeze(alphamatrix(:,:,:,:,mainfold,subfold));
        AUCthreshperfold=squeeze(AUCthreshmatrix(:,:,:,:,mainfold,subfold));
        ROIcritperfold=squeeze(ROIcritmatrix(:,:,:,:,mainfold,subfold));
        [bestvalpersubfold(mainfold, subfold) bestparamspersubfold(mainfold,subfold)]=max(critvarperfold(:));
        findmax=find(critvarperfold==max((critvarperfold(:))));
        
        if length(findmax)==1
            winsub(4,subfold)=alphaperfold(findmax); %these are all the parameter combinations that results in the best outcome...could be one or multiple
            winsub(3,subfold)=lambdaperfold(findmax);
            winsub(2,subfold)=AUCthreshperfold(findmax);
            winsub(1,subfold)=ROIcritperfold(findmax);
        elseif length(findmax)>1
            best_sub=[];
            best_sub(4,:)=alphaperfold(findmax); %these are all the parameter combinations that results in the best outcome...could be one or multiple
            best_sub(3,:)=lambdaperfold(findmax);
            best_sub(2,:)=AUCthreshperfold(findmax);
            best_sub(1,:)=ROIcritperfold(findmax);
            
            winsub(:,subfold)=find_optim(best_sub);
            
        else
            winsub(1,subfold)=NaN;
            winsub(2,subfold)=NaN;
            winsub(3,subfold)=NaN;
            winsub(4,subfold)=NaN;
        end
    end
    if length(find(isnan(winsub(1,:))==1))==length(winsub(1,:));
        disp(['Warning: No parameters found for mainfold ' num2str(mainfold) '. Possibly no variables passed the thresholds in this mainfold.'])
    end
    optim_mainfold(:,mainfold)=find_optim(winsub);
    
    for subfold=1:design.numFolds
        %Stability threshold index
        subfoldidx(1,mainfold,subfold)=optim_mainfold(1,mainfold);
        paramidx_persubfold(1,mainfold,subfold)=winsub(1,subfold);
        if isnan(winsub(1,subfold))==1
            paramidx_persubfold(1,mainfold,subfold)=optim_mainfold(1,mainfold);
        end
        %Prediction error threshold index
        clear diffs
        for zz=1:size(design.merit_thresholds,2)
            diffs(zz)=abs(design.merit_thresholds(mainfold,zz)-winsub(2,subfold));
            diffs2(zz)=abs(design.merit_thresholds(mainfold,zz)-optim_mainfold(2,mainfold));
        end
        ff=find(diffs==min(diffs));
        ff2=find(diffs2==min(diffs2));
        if isnan(winsub(2,subfold))==1
            ff=ff2;
        end
        subfoldidx(2,mainfold,subfold)=ff2(1);
        paramidx_persubfold(2,mainfold,subfold)=ff(1);
        
        %alpha index
        clear diffs
        for zz=1:length(design.alpha)
            diffs(zz)=abs(design.alpha(zz)-winsub(4,subfold));
            diffs2(zz)=abs(design.alpha(zz)-optim_mainfold(4,mainfold));
        end
        ff=find(diffs==min(diffs));
        ff2=find(diffs2==min(diffs2));
        if isnan(winsub(4,subfold))==1
            ff=ff2;
        end
        subfoldidx(4,mainfold,subfold)=ff2(1);
        paramidx_persubfold(4,mainfold,subfold)=ff(1);
        
        %lambda index
        clear diffs
        tmplambdas=squeeze(lambdamatrix(paramidx_persubfold(1,mainfold,subfold),paramidx_persubfold(2,mainfold,subfold),paramidx_persubfold(4,mainfold,subfold),:,mainfold,subfold));
        for zz=1:length(tmplambdas)
            diffs(zz)=abs(tmplambdas(zz)-winsub(3,subfold));
            diffs2(zz)=abs(tmplambdas(zz)-optim_mainfold(3,mainfold));
        end
        ff=find(diffs==min(diffs));
        ff2=find(diffs2==min(diffs2));
        %         try
        paramidx_persubfold(3,mainfold,subfold)=ff(1);
        subfoldidx(3,mainfold,subfold)=ff2(1);
        %         catch ME
        %             paramidx_persubfold(3,mainfold,subfold)=NaN;
        %         end
    end
end
end

function winsub=find_optim(bestsub)
findintersectsub=zeros(size(bestsub));
for j=1:4
    if length(unique(bestsub(j,:)))==length(bestsub(j,:))
        idealval(j)=nanmean(bestsub(j,:));
    elseif length(unique(bestsub(j,:)))==1 %if there is only one value of the parameter
        winsub(j)=bestsub(j,1);
        idealval(j)=bestsub(j,1);
        findintersectsub(j,:)=1;
    else
        idealval(j)=mode(bestsub(j,:)); %mode of all possible parameter values
        if isnan(idealval(j))==1
            idealval(j)=nanmean(bestsub(j,:));
        end
        uq=unique(bestsub(j,:));
        uq(find(isnan(uq)==1))=[];
        clear uqel
        for ll=1:length(uq)
            uqel(ll)=length(find(bestsub(j,:)==uq(ll)));
        end
        tt=find(uqel==max(uqel));
        for ttloop=1:length(tt)
            findintersectsub(j,find(bestsub(j,:)==uq(tt(ttloop))))=1;
        end
    end
end
isctsub=sum(findintersectsub,1);
if length(find(isctsub==max(isctsub)))==1
    winsub=bestsub(:,find(isctsub==max(isctsub)));
else
    go=0;
    f=find(isctsub==max(isctsub));
    while go==0
        for maxloop=1:length(f)
            tmp1=bestsub(:,f(maxloop));
            diff(maxloop)=sum(abs(idealval-tmp1'));
        end
        val2use=find(diff==min(diff));
        try
            winsub=bestsub(:,f(val2use(1)));
            go=1;
        catch ME
            f=find(isctsub==(max(isctsub)-1));
        end
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
