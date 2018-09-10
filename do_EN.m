function [Beta, pred, stats]=do_EN(X, truth, nboot, bagcrit, saveto, type, numFolds, nparam, reg_on)
% reg_on added in July 2017 to see if tuning on specificity rather than
% f1score might make sense with very unbalanced models
% latest update: july 7th 2017
if ~exist(saveto)
    mkdir(saveto)
    save([saveto filesep 'data'], 'X', 'truth', 'nboot')
end
Y=truth;
switch type
    case 'linear'
        family='gaussian';
        link='identity';
        reg_on='x';
    case 'logistic'
        family='binomial';
        link='logit';
        if exist('reg_on', 'var')==0
            %if the groups are balanced
            if length(find(truth==0))/3>length(find(truth==1)) %pos class small --> tune on true pos rate
                reg_on='trueposrate';
            elseif length(find(truth==1))/3>length(find(truth==0)) %neg class small --> tune on true neg rate
                reg_on='truenegrate';
            else
                reg_on='F1score';
            end
        end
end


if exist('nparam', 'var')==0
    nparam=5;
end
lambda= fliplr(logspace(-1.5,0,nparam));
alpha = linspace(0.01,1.0,nparam);

if ~exist([saveto filesep 'LHmerit.mat']) | ~exist([saveto filesep 'lambda_values.mat'])
    [LHmerit lambda_values mf sf]=do_EN_1(X, truth, nboot, bagcrit, saveto, type, nparam,reg_on,numFolds);
else
    load([saveto filesep 'LHmerit.mat'])
    load([saveto filesep 'lambda_values.mat'])
end

for mainfold=1:numFolds
    pause(2)
    
    for subfolds=1:numFolds
        folds=sub2ind([numFolds numFolds], mainfold, subfolds);
        for alpha_looper=1:length(alpha)
            lambda_aggr{mainfold}(subfolds,alpha_looper,:)=squeeze(lambda_values{folds}(alpha_looper, :));
        end
        [r c]=find(squeeze(LHmerit{folds}( :, :))==max(max(squeeze(LHmerit{folds}( :, :)))));
        try
            max_alpha(mainfold, subfolds)=mode(r);
            max_lambda(mainfold, subfolds)=round(median(c(find(r==mode(r)))));
        catch ME
            max_alpha(mainfold, subfolds)=1;
            max_lambda(mainfold, subfolds)=1;
        end
    end
    pause(1)
    alphaidx2use(mainfold)=mode(max_alpha(mainfold,:));
    lambdaidx2use(mainfold)=mode(max_lambda(mainfold,:));
    alpha2use(mainfold)=alpha(alphaidx2use(mainfold));
    if lambdaidx2use(mainfold)<nparam
        lambda2use(mainfold)=mean(squeeze(lambda_aggr{mainfold}(:,alphaidx2use(mainfold),lambdaidx2use(mainfold)+1)));
    else
        lambda2use(mainfold)=mean(squeeze(lambda_aggr{mainfold}(:,alphaidx2use(mainfold),lambdaidx2use(mainfold))));
    end
    
    trainsubjects = find(mf ~= mainfold);
    testsubjects=find(mf == mainfold);
    
    options=glmnetSet;
    options.standardize=true;
    options.alpha=alpha2use(mainfold);
    options.lambda=lambda2use(mainfold);
    
    fit=glmnet(X(trainsubjects,:),truth(trainsubjects),family,options);
    B0=fit.beta;
    Beta{mainfold}=[fit.a0; B0];
    GetProbs{mainfold} = glmval(Beta{mainfold},X(testsubjects,:),link);
    pred(testsubjects)=GetProbs{mainfold};
    
    %% to test overfit do pred for trainsubs also
    tstGetProbs{mainfold} = glmval(Beta{mainfold},X(trainsubjects,:),link);
    tstpred(trainsubjects)=tstGetProbs{mainfold};
    %     switch type
    %         case 'linear',
    %             [Merit.train_r(mainfold),Merit.train_p(mainfold)]=corr(tstpred(trainsubjects)', truth(trainsubjects));
    %             Merit.train_mse(mainfold)=(abs(truth(trainsubjects)-tstpred(trainsubjects)')'*abs(truth(trainsubjects)-tstpred(trainsubjects)')/length(truth(trainsubjects)));
    %         case 'logistic'
    %             %             [Merit.train_overall_AUC{mainfold},Merit.train_fpr{mainfold},Merit.train_tpr{mainfold}] = fastAUC(logical(truth),tstpred',0);
    %             %             [Merit.train_prec{mainfold}, Merit.train_tpr{mainfold}, Merit.train_fpr{mainfold}, Merit.train_thresh{mainfold}] = prec_rec_rob_mod(tstpred', logical(truth),'tstPrecRec', 'plotPR',0);
    %             %             fscoretrain_=(Merit.train_prec{mainfold}.*Merit.train_tpr{mainfold})./(Merit.train_prec{mainfold}+Merit.train_tpr{mainfold});
    %             %             Merit.train_F1score{mainfold}=max(fscoretrain_);
    %             try
    %                 [Merit.train_accuracy{mainfold} Merit.train_sensitivity{mainfold} Merit.train_specificity{mainfold} Merit.train_AUC{mainfold} Merit.train_F1{mainfold}]=class_mets_new(logical(truth), tstpred);
    %             catch
    %                 disp('oh well')
    %             end
    %     end
end

switch type
    case 'linear',
        [stats.r, stats.p]=corr(pred', truth);
        stats.mse=(abs(truth-pred')'*abs(truth-pred')/length(truth));
    case 'logistic'
        %         [stats.overall_AUC,stats.fpr,stats.tpr] = fastAUC(logical(truth),pred',0);
        %         [stats.prec, stats.tpr, stats.fpr, stats.thresh] = prec_rec_rob_mod(pred', logical(truth),'tstPrecRec', 'plotPR',0);
        %         fscore=(stats.prec.*stats.tpr)./(stats.prec+stats.tpr);
        %         stats.F1score=max(fscore);
        [stats.accuracy stats.trueposrate stats.truenegrate stats.falseposrate stats.falsenegrate stats.precision stats.NPV stats.AUC stats.F1score]=class_mets_new(truth, pred);
end


if ~exist(saveto, 'dir')
    mkdir(saveto);
end
cd(saveto)
save('data', 'X', 'truth');
save('betas', 'Beta');
save('prediction', 'pred');
save('params2use', 'alpha2use', 'lambda2use');
save('results', 'stats');
% save('Merit',  'Merit');
save('Results',  'pred', 'Beta', 'alpha2use', 'lambda2use', 'stats', 'mf', 'sf');
end

function [LHmerit lambda_values mf sf]=do_EN_1(X, truth, nboot, bagcrit, saveto, type, nparam, reg_on, numFolds)

% latest update: may 30th 2017

Y=truth;
switch type
    case 'linear'
        family='gaussian';
        link='identity';
    case 'logistic'
        family='binomial';
        link='logit';
end
ok=0;
while ok==0
    [mf, sf]=AssignFolds(size(X,1),numFolds,numFolds);
    [mf ,sf]=check_CV(mf, sf,truth, numFolds);
    folds=0;
    while folds<numFolds*numFolds
        folds=folds+1;
        tmplambda2use=[];
        [mainfold, subfolds]=ind2sub([numFolds numFolds], folds);
        trainsubjectsTune = find(mf ~= mainfold & sf(:,mainfold) ~=subfolds);
        testsubjectsTune = find(mf ~= mainfold & sf(:,mainfold)==subfolds);
        if length(unique(truth(testsubjectsTune)))<2 | length(unique(truth(trainsubjectsTune)))<2
            folds=numFolds*numFolds;
            ok=0;
        elseif folds==numFolds*numFolds
            ok=1;
        end
    end
end


lambda= fliplr(logspace(-1.5,0,nparam));
alpha = linspace(0.01,1.0,nparam);

%%PARFOR
parfor folds=1:numFolds*numFolds
%     disp(folds)
    tmplambda2use=[];
    [mainfold, subfolds]=ind2sub([numFolds numFolds], folds);
    trainsubjectsTune = find(mf ~= mainfold & sf(:,mainfold) ~=subfolds);
    testsubjectsTune = find(mf ~= mainfold & sf(:,mainfold)==subfolds);
    
    options=glmnetSet;
    options.standardize=true;
    options.nlambda=length(lambda);
    options.lambda_min=min(lambda);
    for alpha_looper=1:length(alpha)
        options.alpha=alpha(alpha_looper);
        tmpmerit=[];
        tmpvars2use=[];
        if nboot<2
            fit=glmnet(X(trainsubjectsTune, :),truth(trainsubjectsTune),family,options);
            B0=fit.beta; intercept=fit.a0;
            while size(B0,2)<options.nlambda
                options.lambda=linspace(lambda(1)+lambda(1)/5, 1, length(lambda));
                fit=glmnet(X(trainsubjectsTune, :),truth(trainsubjectsTune),family,options);
                B0=fit.beta; intercept=fit.a0;
            end
            lambda_values{folds}(alpha_looper, :)=fit.lambda;
            if size(intercept,2)==size(B0,2)
                b=[intercept; B0];
            elseif size(intercept,1)==size(B0,2)
                b=[intercept'; B0];
            elseif size(intercept,2)==size(B0,1)
                b=[intercept; B0'];
            elseif size(intercept,1)==size(B0,1)
                b=[intercept'; B0'];
            else
                disp(['could not concatenate intercept and B0.'])
                size(intercept)
                size(B0)
            end
            
            for lambda_looper=1:length(lambda)
                bb=b(:,lambda_looper);
                getprobstune = glmval(squeeze(bb),X(testsubjectsTune,:),link);
                switch type
                    case 'linear',
                        LHmerit{folds} (alpha_looper,lambda_looper) = -sqrt(abs(truth(testsubjectsTune)-getprobstune)'*abs(truth(testsubjectsTune)-getprobstune)/length(truth(testsubjectsTune)));
                    case 'logistic',
%                         try
                            [accuracy trueposrate truenegrate falseposrate falsenegrate precision NPV AUC F1score]=class_mets_new(truth(testsubjectsTune), getprobstune');
%                         catch
%                             pause
%                         end
                        if strcmp(reg_on, 'trueposrate')==1
                            LHmerit{folds} (alpha_looper,lambda_looper) =trueposrate;
                        elseif strcmp(reg_on, 'truenegrate')==1
                            LHmerit{folds} (alpha_looper,lambda_looper) =truenegrate;
                        elseif strcmp(reg_on, 'accuracy')==1
                            LHmerit{folds} (alpha_looper,lambda_looper) =accuracy;
                        elseif strcmp(reg_on, 'AUC')==1
                            LHmerit{folds} (alpha_looper,lambda_looper) =AUC;
                        elseif strcmp(reg_on, 'F1score')==1
                            LHmerit{folds} (alpha_looper,lambda_looper) =F1score;
                        elseif strcmp(reg_on, 'precision')==1
                            LHmerit{folds} (alpha_looper,lambda_looper) =posretrieval;
                        end
                end
            end
        elseif nboot>1
            for bootct=1:nboot
                train=0;
                while length(train)<2
                    [Xboot,Yboot,indexselect]=bootstrapal(X(trainsubjectsTune,:),truth(trainsubjectsTune),2/3);
                    train=unique(Yboot);
                end
                
                fit=glmnet(Xboot,Yboot,family,options);
                B0=fit.beta; intercept=fit.a0;
                
                while size(B0,2)<options.nlambda
                    options.lambda=linspace(lambda(1)+lambda(1)/10, 1, length(lambda));
                    fit=glmnet(Xboot,Yboot,family,options);
                    B0=fit.beta; intercept=fit.a0;
                end
                try
                    b=[intercept'; B0];
                catch
                    try
                        b=[intercept; B0];
                    catch
                        try
                            b=[intercept'; B0'];
                        catch
                            b=[intercept;B0'];
                        end
                    end
                end
                tmplambda2use(bootct,:)=fit.lambda;
                for lambda_looper=1:length(lambda)
                    bb=b(:,lambda_looper);
                    getprobstune = glmval(squeeze(bb),X(testsubjectsTune,:),link);
                    switch type
                        case 'linear',
                            tmp1 = -sqrt(abs(Y(testsubjectsTune)-getprobstune)'*abs(Y(testsubjectsTune)-getprobstune)/length(Y(testsubjectsTune)));
                        case 'logistic'
                            try
                                [accuracy trueposrate truenegrate falseposrate falsenegrate posretrieval negretrieval AUC F1score]=class_mets_new(Y(testsubjectsTune), getprobstune');
                            catch
                                pause
                            end
                            if strcmp(reg_on, 'trueposrate')==1
                                tmp1 =trueposrate;
                            elseif strcmp(reg_on, 'truenegrate')==1
                                tmp1 =truenegrate;
                            elseif strcmp(reg_on, 'accuracy')==1
                                tmp1 =accuracy;
                            elseif strcmp(reg_on, 'AUC')==1
                                tmp1 =AUC;
                            elseif strcmp(reg_on, 'F1score')==1
                                tmp1 =F1score;
                            elseif strcmp(reg_on, 'posretrieval')==1
                                tmp1 =posretrieval;
                            end
                    end
                    tmpmerit(bootct,lambda_looper)=tmp1;
                end
            end
            lambda_values{folds}(alpha_looper, :)=mean(tmplambda2use);
            for lambda_looper=1:length(lambda)
                if strcmp('cdf', bagcrit)==1
                    tmpy=cdf('norm', tmpmerit(:,lambda_looper), mean(tmpmerit(:,lambda_looper)), std(tmpmerit(:,lambda_looper)));
                    findLH=tmpy(find(tmpy<=0.05));
                    findLH=max(findLH);
                    if isempty(findLH)==0
                        tmpLHmerit=tmpmerit(find(tmpy==findLH),lambda_looper);
                        tmpLH_y=findLH;
                    else
                        tmpLHmerit=NaN;
                        tmpLH_y=NaN;
                    end
                    LHmerit{folds}(alpha_looper, lambda_looper)=tmpLHmerit(1);
                elseif strcmp('median', bagcrit)==1
                    LHmerit{folds}(alpha_looper, lambda_looper)=median(tmpmerit(:,lambda_looper));
                else
                    disp('please enter a valid bootstrap aggregation criterion (cdf or mean)');
                end
            end
        end
    end
end
save([saveto filesep 'LHmerit.mat'], 'LHmerit', 'mf', 'sf')
save([saveto filesep 'lambda_values.mat'], 'lambda_values')
end
