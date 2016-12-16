function [merit_per_var] = RAFT_FS(design)
% feature selection step of Regularized Adaptive Feature Thresholding
%
% for comments and questions please contact lee.jollans@gmail.com

cd(design.saveto);

fprintf('Performing feature thresholding with %d subjects and %d predictors\n', size(design.data));

design.nvars=size(design.data,2);
if design.nboot>1
    for n=1:design.nboot
        tmpmerit{n}=NaN(design.numFolds*design.numFolds,design.nvars);
    end
    for bootct=1:design.nboot
        for folds=1:(design.numFolds*design.numFolds)
            tmp_merit{folds}=NaN(design.nvars,1);
        end
        parfor folds=1:(design.numFolds*design.numFolds)
            [outerFold, middleFold]=ind2sub([design.numFolds design.numFolds], folds);
            trainingsubs=find(design.subfolds(:,outerFold)~=middleFold & design.subfolds(:,outerFold)~=-1); 
            testsubs=find(design.subfolds(:,outerFold)==middleFold & design.subfolds(:,outerFold)~=-1);
            try
                [Xboot,Yboot]=bootstrapal(design.data([trainingsubs; testsubs],:),design.outcome([trainingsubs; testsubs]),design.Ratio);
            catch ME
                [Xboot,Yboot]=bootstrapal(design.data([trainingsubs; testsubs],:),design.outcome([trainingsubs; testsubs])',design.Ratio);
            end
            for vars=1:design.nvars
                
                switch(design.type)
                    case 'linear',
                        [b,dev,stats]=glmfit(Xboot(1:length(trainingsubs),vars),Yboot(1:length(trainingsubs)),design.distribution,'link',design.link);
                        pred=glmval(b,[Xboot(length(trainingsubs)+1:end,vars)],design.link, 'constant', 'on');
                        truth=Yboot(length(trainingsubs)+1:end);
                        tmp = -sqrt(abs(truth-pred)'*abs(truth-pred)/length(truth));
                        tmp_merit{folds}(vars)=tmp;
                    case 'logistic',
                        [b,dev,stats]=glmfit(Xboot(1:length(trainingsubs),vars),Yboot(1:length(trainingsubs)),design.distribution,'link',design.link);
                        pred = glmval(b,Xboot(length(trainingsubs)+1:end,vars)',design.link, 'constant', 'on');
                        truth=Yboot(length(trainingsubs)+1:end);
                        switch(design.balanced)
                            case 'balanced'
                                [tmp1,fpr,tpr] = fastAUC(truth,pred,0);
                                tmp_merit{folds}(vars)=tmp1;
                            case 'unbalanced'
                                [prec, tpr, fpr, thresh] = prec_rec_rob_mod(pred, truth,'tstPrecRec', 'plotPR',0, 'numThresh',100);
                                fscore=(prec.*tpr)./(prec+tpr);
                                tmp_merit{folds}(vars)=max(fscore);
                        end
                end
            end
        end
        for folds=1:(design.numFolds*design.numFolds)
            tmpmerit{bootct}(folds,:)=tmp_merit{folds};
        end
    end
    for n=1:design.numFolds*design.numFolds
        merit_per_var{n}=NaN(design.nvars,1);
    end
    parfor folds=1:design.numFolds*design.numFolds
        for vars=1:design.nvars
            x=[];
            for bootct=1:design.nboot
                x=[x,tmpmerit{bootct}(folds, vars)];
            end
            if strcmp(design.bagcrit, 'cdf')==1
                tmpy=cdf('norm', x, mean(x), std(x));
                findLH=tmpy(find(tmpy<=design.siglevel));
                findLH=max(findLH);
                if isempty(findLH)==0
                    tmpLHmerit=x(find(tmpy==findLH));
                    tmpLH_y=findLH;
                else
                    tmpLHmerit=NaN;
                    tmpLH_y=NaN;
                end
                merit_per_var{folds}(vars)=tmpLHmerit(1);
            elseif strcmp(design.bagcrit, 'median')==1
                merit_per_var{folds}(vars)=median(x);
            else
                disp('Please enter a valid method for bagging (mean or cdf)')
            end
        end
    end
else
    for n=1:(design.numFolds*design.numFolds)
        merit_per_var{n}=NaN(design.nvars,1);
    end
    parfor folds=1:(design.numFolds*design.numFolds)
        [outerFold, middleFold]=ind2sub([design.numFolds design.numFolds], folds);
        trainingsubs=find(design.subfolds(:,outerFold)~=middleFold & design.subfolds(:,outerFold)~=-1); %subs in the training set for each inner fold
        testsubs=find(design.subfolds(:,outerFold)==middleFold & design.subfolds(:,outerFold)~=-1);
        for vars=1:design.nvars
            
            switch(design.type)
                case 'linear',
                    [b,dev,stats]=glmfit(design.data(trainingsubs,vars),design.outcome(trainingsubs),design.distribution,'link',design.link);
                    pred=glmval(b,[design.data(testsubs,vars)],design.link, 'constant', 'on');
                    truth=design.outcome(testsubs);
                    tmp = -sqrt(abs(truth-pred)'*abs(truth-pred)/length(truth));
                    merit_per_var{folds}(vars)=tmp;
                case 'logistic',
                    [b,dev,stats]=glmfit(design.data(trainingsubs,vars),design.outcome(trainingsubs),design.distribution,'link',design.link);
                    pred = glmval(b,design.data(testsubs,vars)',design.link, 'constant', 'on');
                    truth=design.outcome(testsubs);
                    switch(design.balanced)
                        case 'balanced'
                            [tmp1,fpr,tpr] = fastAUC(truth,pred,0);
                            merit_per_var{folds}(vars)=tmp1;
                        case 'unbalanced'
                            [prec, tpr, fpr, thresh] = prec_rec_rob_mod(pred, truth,'tstPrecRec', 'plotPR',0, 'numThresh',100);
                            fscore=(prec.*tpr)./(prec+tpr);
                            merit_per_var{folds}(vars)=max(fscore);
                    end
            end
        end
    end
end

cd(design.saveto);
save('merit_per_var', 'merit_per_var');
end
