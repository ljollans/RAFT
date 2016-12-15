function RAFT_2nd_Level(design, pass_vars, folds2run)
% Second level for Regularized Adaptive Feature Thresholding
% Here models for each combination of parameters and thresholds are
% evaluated and the superior models are selected at the third level
%
% for comments and questions please contact lee.jollans@gmail.com

cd(design.saveto);

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
            ttmpvars2use=zeros(design.nboot,design.numLambdas,design.nvars);
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
                    while length(train)<2
                        %get bootstrap data subset
                        try
                        [Xboot,Yboot]=bootstrapal([design.data([trainsubjectsTune; testsubjectsTune],:), design.extradata([trainsubjectsTune; testsubjectsTune],:)],design.outcome([trainsubjectsTune; testsubjectsTune]),design.Ratio);
                        catch ME
                            [Xboot,Yboot]=bootstrapal([design.data([trainsubjectsTune; testsubjectsTune],:), design.extradata([trainsubjectsTune; testsubjectsTune],:)],design.outcome([trainsubjectsTune; testsubjectsTune])',design.Ratio);
                        end
                        %make sure we have multiple training subjects
                        train=unique(Yboot(length(trainsubjectsTune)+1:end));
                    end
                    
                    %fit the elastic net and get beta and intercept values
                    fit=glmnet(Xboot(1:length(trainsubjectsTune), Vars2pick),Yboot(1:length(trainsubjectsTune)),design.family,options);
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
                    usedvars=zeros(design.numLambdas,design.nvars);
                    
                    for lambda_looper=1:length(design.lambda)
                        bb=b(:,lambda_looper);
                        %identify what variables the elastic net excluded
                        tmp=zeros(design.nvars,1);
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
                    usedvars=zeros(design.numLambdas,design.nvars);
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
    save([design.saveto filesep 'Second_level_results_' num2str(mainfold) '_' num2str(subfolds)], 'LHmerit2save', 'vars2use2save', 'lambda_values2save');
end
