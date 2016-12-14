function [params2pick, Merit, Beta, Vars2pick_main,  GetProbs, design, stats] = RAFT_Model_Selection(design, vars2use, LHmerit, lambda_values, pass_vars, use_pass_vars, startROIthresh, startmeritthresh)
% select the best performing model across subfolds for each main CV fold
% and applies it to the test set
%
% for comments and questions please contact lee.jollans@gmail.com

design.prediction=[];
fprintf('Performing Mainfold Analysis with %d subjects and %d predictors\n', size(design.data));

outcome2use=design.outcome;

options=glmnetSet;
options.standardize=true;
options.nlambda=1;

[params2pick paramidx lambdavalues] = findbestparams(LHmerit, design, lambda_values, startROIthresh, startmeritthresh);
intalpha=4; intlambda=3; intmerit=2; introi=1;
params2use=params2pick;
paramidx2use=paramidx;
for mainfold=1:design.numFolds
    %%
    Vars2pick_main{mainfold}=[];
    chkr=0;
    
    while length(Vars2pick_main{mainfold})<1
        paramidx_perfold=squeeze(paramidx2use(:,mainfold,:));
        [rf cf]=find(isnan(paramidx_perfold)==1);
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
        
        if length(Vars2pick_main{mainfold})<1
            for roiidxs=1:design.numFolds
                for meritidxs=1:size(design.merit_thresholds,2)
                    for alphaidxs=1:design.numAlphas
                        for lambdaidxs=1:design.numLambdas
                            for subfold=1:design.numFolds
                                folds=sub2ind([design.numFolds design.numFolds], mainfold, subfold);
                                tmptrackvars(subfold,:)=squeeze(vars2use(mainfold,subfold,roiidxs, meritidxs, alphaidxs, lambdaidxs,:));
                            end
                            varscount=sum(tmptrackvars);
                            pickr(roiidxs, meritidxs, alphaidxs, lambdaidxs)=not(isempty(find(varscount>=roiidxs)));
                        end
                    end
                end
            end
            while length(Vars2pick_main{mainfold})<1
                pickr(tmpidx(introi),tmpidx(intmerit),tmpidx(intalpha),tmpidx(intlambda))=0;
                
                if sum(sum(sum(sum(pickr))))==0
                    use_pass_vars=1;
                    cd(design.saveto)
                    load('pass_vars')
                    for roiidxs=1:design.numFolds
                        for meritidxs=1:size(design.merit_thresholds,2)
                            for alphaidxs=1:design.numAlphas
                                for lambdaidxs=1:design.numLambdas
                                    for subfold=1:design.numFolds
                                        tmptrackvars(subfold,:)=squeeze(pass_vars(roiidxs, meritidxs, mainfold, subfold,:));
                                    end
                                    varscount=sum(tmptrackvars);
                                    pickr(roiidxs, meritidxs, alphaidxs, lambdaidxs)=not(isempty(find(varscount>=roiidxs)));
                                end
                            end
                        end
                    end
                end
                
                [roiidx meritidx alphaidx lambdaidx]=ind2sub([design.numFolds size(design.merit_thresholds,2) design.numAlphas design.numLambdas],find(pickr));
                tmpoptimalchk=zeros(design.numFolds,length(roiidx),4);
                for subfold=1:design.numFolds
                    tmpoptimalchk(subfold,find(roiidx==paramidx2use(introi,mainfold,subfold)), introi)=1;
                    tmpoptimalchk(subfold,find(meritidx==paramidx2use(intmerit,mainfold,subfold)), intmerit)=1;
                    tmpoptimalchk(subfold,find(lambdaidx==paramidx2use(intlambda,mainfold,subfold)), intlambda)=1;
                    tmpoptimalchk(subfold,find(alphaidx==paramidx2use(intalpha,mainfold,subfold)), intalpha)=1;
                end
                optimalchk=squeeze(sum(sum(tmpoptimalchk,3)));
                tt=find(optimalchk==max(optimalchk));
                clear diff_params tmpdiff_params
                for diffloop=1:length(tt)
                    for subfold=1:design.numFolds
                        tmpdiff_params(diffloop,subfold)=sum(abs([roiidx(tt(diffloop)) meritidx(tt(diffloop)) lambdaidx(tt(diffloop)) alphaidx(tt(diffloop))]-paramidx_perfold(:,subfold)'));
                    end
                    diff_params(diffloop)=sum(tmpdiff_params(diffloop,:));
                end
                paramidxs=tt(find(diff_params==min(diff_params)));
                for subfold=1:design.numFolds
                    paramidx2use(introi,mainfold,subfold)=roiidx(paramidxs(1));
                    paramidx2use(intmerit,mainfold,subfold)=meritidx(paramidxs(1));
                    paramidx2use(intalpha,mainfold,subfold)=alphaidx(paramidxs(1));
                    paramidx2use(intlambda,mainfold,subfold)=lambdaidx(paramidxs(1));
                    lambdacollect(subfold)=lambdavalues(roiidx(paramidxs(1)),meritidx(paramidxs(1)),alphaidx(paramidxs(1)),lambdaidx(paramidxs(1)),mainfold,subfold);
                end
                params2use(introi,mainfold)=roiidx(paramidxs(1));
                params2use(intmerit,mainfold)=design.merit_thresholds(mainfold,meritidx(paramidxs(1)));
                params2use(intalpha,mainfold)=design.alpha(alphaidx(paramidxs(1)));
                params2use(intlambda,mainfold)=mean(lambdacollect);
                
                paramidx_perfold=squeeze(paramidx2use(:,mainfold,:));
                [rf cf]=find(isnan(paramidx_perfold)==1);
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
                        [prec, tpr, fpr, thresh] = prec_rec_rob_mod_ej(GetProbs{mainfold}, logical(design.outcome(testsubjects)),'tstPrecRec', 'plotPR',0);
                        [Merit.AUC_unbalanced(mainfold),fpr2,tpr2] = fastAUC(logical(design.outcome(testsubjects)),GetProbs{mainfold},0);
                    catch ME
                        length(testsubjects);
                    end
                    fscore=(prec.*tpr)./(prec+tpr);
                    Merit.F1score(mainfold)=max(fscore);
            end
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
        [stats.prec, stats.tpr, stats.fpr, stats.thresh] = prec_rec_rob_mod_ej(design.prediction', logical(design.outcome),'tstPrecRec', 'plotPR',0);
        fscore=(stats.prec.*stats.tpr)./(stats.prec+stats.tpr);
        stats.F1score=max(fscore);
        save('Results',  'Merit', 'GetProbs', 'Beta', 'params2pick', 'params2use', 'stats', 'Vars2pick_main', 'design');
end
end

function [optim_mainfold paramidx_persubfold lambdamatrix] = findbestparams(AUCTune_in, design, lambda_values, startROIthresh, startAUCthresh)

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
        %Stability trheshold index
        paramidx_persubfold(1,mainfold,subfold)=winsub(1,subfold);
        if isnan(winsub(1,subfold))==1
            paramidx_persubfold(1,mainfold,subfold)=optim_mainfold(1,mainfold);
        end
        %Prediction error threshold index
        clear diffs
        for zz=1:size(design.merit_thresholds,2)
            diffs(zz)=abs(design.merit_thresholds(mainfold,zz)-winsub(2,subfold));
        end
        ff=find(diffs==min(diffs));
        if isnan(winsub(2,subfold))==1
            for zz=1:size(design.merit_thresholds,2)
            diffs(zz)=abs(design.merit_thresholds(mainfold,zz)-optim_mainfold(2,mainfold));
        end
        ff=find(diffs==min(diffs));
            paramidx_persubfold(2,mainfold,subfold)=ff(1);
        else
        paramidx_persubfold(2,mainfold,subfold)=ff(1);
        end

        %alpha index
        clear diffs
        for zz=1:length(design.alpha)
            diffs(zz)=abs(design.alpha(zz)-winsub(4,subfold));
        end
        ff=find(diffs==min(diffs));
        if isnan(winsub(4,subfold))==1
                    for zz=1:length(design.alpha)
            diffs(zz)=abs(design.alpha(zz)-optim_mainfold(4,mainfold));
        end
        ff=find(diffs==min(diffs));
            paramidx_persubfold(4,mainfold,subfold)=ff(1);
        else
        paramidx_persubfold(4,mainfold,subfold)=ff(1);
        end
        %lambda index
        clear diffs
        tmplambdas=squeeze(lambdamatrix(paramidx_persubfold(1,mainfold,subfold),paramidx_persubfold(2,mainfold,subfold),paramidx_persubfold(4,mainfold,subfold),:,mainfold,subfold));
        for zz=1:length(tmplambdas)
            diffs(zz)=abs(tmplambdas(zz)-winsub(3,subfold));
        end
        ff=find(diffs==min(diffs));
        try
        paramidx_persubfold(3,mainfold,subfold)=ff(1);
        catch ME
            paramidx_persubfold(3,mainfold,subfold)=NaN;
        end
    end
end
end

function winsub=find_optim(bestsub)
findintersectsub=zeros(size(bestsub));
for j=1:4
    if length(unique(bestsub(j,:)))==length(bestsub(j,:))
        idealval(j)=nanmean(bestsub(j,:));
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
    if length(unique(bestsub(j,:)))==1 %if there is only one value of the parameter
        winsub(j)=bestsub(j,1);
        findintersectsub(j,:)=1;
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