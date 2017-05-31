function [Beta, pred, stats]=do_EN(X, truth, nboot, bagcrit, saveto, type)

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
[mf, sf]=AssignFolds(size(X,1),10,10);

lambda= fliplr(logspace(-3,0,5));
alpha = linspace(0.01,1.0,5);

for folds=1:100
    tmplambda2use=[];
    [mainfold subfolds]=ind2sub([10 10], folds);
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
            lambda_values{folds}(alpha_looper, :)=fit.lambda;
            while size(B0,2)<options.nlambda
                options.lambda=linspace(lambda(1)+lambda(1)/10, 1, length(lambda));
                fit=glmnet(X(trainsubjectsTune, :),truth(trainsubjectsTune),family,options);
                B0=fit.beta; intercept=fit.a0;
                lambda_values{folds}(alpha_looper, :)=fit.lambda;
            end
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
                        [prec, tpr] = prec_rec_rob_mod(getprobstune, truth(testsubjectsTune),'tstPrecRec', 'plotPR',0, 'numThresh',100);
                        fscore=(prec.*tpr)./(prec+tpr);
                        LHmerit{folds} (alpha_looper,lambda_looper) =max(fscore);
                end
            end
        elseif nboot>1
            for bootct=1:nboot
                train=0;
                while length(train)<2
                    [Xboot,Yboot,indexselect]=bootstrapal(X(trainsubjectsTune,:),truth,2/3);
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
                        case 'logistic',
                            [prec, tpr] = prec_rec_rob_mod(getprobstune, Y(testsubjectsTune),'tstPrecRec', 'plotPR',0, 'numThresh',100);
                            fscore=(prec.*tpr)./(prec+tpr);
                            tmp1=max(fscore);
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


for mainfold=1:10
    for subfolds=1:10
        folds=sub2ind([10 10], mainfold, subfolds);
        for alpha_looper=1:length(alpha)
            lambda_aggr{mainfold}(subfolds,alpha_looper,:)=squeeze(lambda_values{folds}(alpha_looper, :));
        end
        [r c]=find(squeeze(LHmerit{folds}( :, :))==max(max(squeeze(LHmerit{folds}( :, :)))));
        try
            max_alpha(mainfold, subfolds)=r(1);
            max_lambda(mainfold, subfolds)=c(1);
        catch ME
            max_alpha(mainfold, subfolds)=1;
            max_lambda(mainfold, subfolds)=1;
        end
    end
    alphaidx2use(mainfold)=mode(max_alpha(mainfold,:));
    lambdaidx2use(mainfold)=mode(max_lambda(mainfold,:));
    alpha2use(mainfold)=alpha(alphaidx2use(mainfold));
    lambda2use(mainfold)=mean(squeeze(lambda_aggr{mainfold}(:,alphaidx2use(mainfold),lambdaidx2use(mainfold))));
    
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
    switch type
        case 'linear',
            [Merit.train_r(mainfold),Merit.train_p(mainfold)]=corr(tstpred(trainsubjects)', truth(trainsubjects));
            Merit.train_mse(mainfold)=(abs(truth(trainsubjects)-tstpred(trainsubjects)')'*abs(truth(trainsubjects)-tstpred(trainsubjects)')/length(truth(trainsubjects)));
        case 'logistic'
            [Merit.overall_AUC{mainfold},Merit.fpr{mainfold},Merit.tpr{mainfold}] = fastAUC(logical(truth),pred',0);
            [Merit.prec{mainfold}, Merit.tpr{mainfold}, Merit.fpr{mainfold}, Merit.thresh{mainfold}] = prec_rec_rob_mod(pred', logical(truth),'tstPrecRec', 'plotPR',0);
            fscore=(Merit.prec{mainfold}.*Merit.tpr{mainfold})./(Merit.prec{mainfold}+Merit.tpr{mainfold});
            Merit.F1score{mainfold}=max(fscore);
    end
end

switch type
    case 'linear',
        [stats.r, stats.p]=corr(pred', truth);
        stats.mse=(abs(truth-pred')'*abs(truth-pred')/length(truth));
    case 'logistic'
        [stats.overall_AUC,stats.fpr,stats.tpr] = fastAUC(logical(truth),pred',0);
        [stats.prec, stats.tpr, stats.fpr, stats.thresh] = prec_rec_rob_mod(pred', logical(truth),'tstPrecRec', 'plotPR',0);
        fscore=(stats.prec.*stats.tpr)./(stats.prec+stats.tpr);
        stats.F1score=max(fscore);
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
% save('Merit', , 'Merit');
save('Results',  'pred', 'Beta', 'alpha2use', 'lambda2use', 'stats', 'mf', 'sf');
end
