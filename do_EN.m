function [Beta, pred, r, p, mse]=do_EN(X, truth, nboot, bagcrit, saveto)
[mf, sf]=AssignFolds(size(X,1),10,10);

lambda= fliplr(logspace(-3,0,5));
alpha = linspace(0.01,1.0,5);

parfor folds=1:100
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
            fit=glmnet(X(trainsubjectsTune, :),truth(trainsubjectsTune),'gaussian',options);
            B0=fit.beta; intercept=fit.a0;
            lambda_values{folds}(alpha_looper, :)=fit.lambda;
            while size(B0,2)<options.nlambda
                options.lambda=linspace(lambda(1)+lambda(1)/10, 1, length(lambda));
                fit=glmnet(X(trainsubjectsTune, :),truth(trainsubjectsTune),'gaussian',options);
                B0=fit.beta; intercept=fit.a0;
                lambda_values{folds}(alpha_looper, :)=fit.lambda;
            end
            b=[intercept'; B0];
            for lambda_looper=1:length(lambda)
                bb=b(:,lambda_looper);
                getprobstune = glmval(squeeze(bb),X(testsubjectsTune,:),'identity');
                LHmerit{folds} (alpha_looper,lambda_looper) = -sqrt(abs(truth(testsubjectsTune)-getprobstune)'*abs(truth(testsubjectsTune)-getprobstune)/length(truth(testsubjectsTune)));
            end
        elseif nboot>1
            for bootct=1:nboot
                train=0;
                while length(train)<2
                    [Xboot,Yboot,indexselect]=bootstrapal(X,truth,2/3);
                    train=unique(Yboot(testsubjectsTune));
                end
                
                fit=glmnet(Xboot(trainsubjectsTune, :),Yboot(trainsubjectsTune),'gaussian',options);
                B0=fit.beta; intercept=fit.a0;
                
                while size(B0,2)<options.nlambda
                    options.lambda=linspace(lambda(1)+lambda(1)/10, 1, length(lambda));
                    fit=glmnet(Xboot(trainsubjectsTune, :),Yboot(trainsubjectsTune),'gaussian',options);
                    B0=fit.beta; intercept=fit.a0;
                end
                b=[intercept'; B0];
                tmplambda2use(bootct,:)=fit.lambda;
                for lambda_looper=1:length(lambda)
                    bb=b(:,lambda_looper);
                    getprobstune = glmval(squeeze(bb),Xboot(testsubjectsTune,:),'identity');
                    tmp1 = -sqrt(abs(Yboot(testsubjectsTune)-getprobstune)'*abs(Yboot(testsubjectsTune)-getprobstune)/length(Yboot(testsubjectsTune)));
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
    
    options.alpha=alpha2use(mainfold);
    options.lambda=lambda2use(mainfold);
    
    fit=glmnet(X(trainsubjects,:),truth(trainsubjects),'gaussian',options);
    B0=fit.beta;
    Beta{mainfold}=[fit.a0; B0];
    GetProbs{mainfold} = glmval(Beta{mainfold},X(testsubjects,:),'identity');
    pred(testsubjects)=GetProbs{mainfold};
end
[r, p]=corr(pred', truth);
mse=(abs(truth-pred')'*abs(truth-pred')/length(truth));
   
if ~exist(saveto, 'dir')
    mkdir(saveto);
end
cd(saveto)
save('data', 'X', 'truth');
save('betas', 'Beta');
save('prediction', 'pred');
save('params2use', 'alpha2use', 'lambda2use');
save('results', 'r', 'p', 'mse');
save('Results',  'pred', 'Beta', 'alpha2use', 'lambda2use', 'r', 'p', 'mse');
end