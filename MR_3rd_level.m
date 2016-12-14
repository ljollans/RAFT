function [params2pick, Merit, Beta, Vars2pick_main,  GetProbs, design, stats] = MR_3rd_level(design, vars2use, LHmerit, startROIthresh, startmeritthresh)

design.prediction=[];
fprintf('Performing Mainfold Analysis with %d subjects and %d predictors\n', size(design.data));

outcome2use=design.outcome;

[params2pick paramidx] = findbestparams_MR(LHmerit, design, startROIthresh, startmeritthresh);
intmerit=2; introi=1;
params2use=params2pick;
paramidx2use=paramidx;
for mainfold=1:design.numFolds
    use_pass_vars=0;
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
                    trackvars(mainfold,subfold,:)=squeeze(vars2use(mainfold, subfold,tmpidx(introi),tmpidx(intmerit),:));
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
                for meritidxs=1:length(design.merit_thresholds)
                            for subfold=1:design.numFolds
                                folds=sub2ind([design.numFolds design.numFolds], mainfold, subfold);
                                tmptrackvars(subfold,:)=squeeze(vars2use(mainfold,subfold,roiidxs, meritidxs,:));
                            end
                            varscount=sum(tmptrackvars);
                            pickr(roiidxs, meritidxs)=not(isempty(find(varscount>=roiidxs)));
                end
            end
            while length(Vars2pick_main{mainfold})<1
                pickr(tmpidx(introi),tmpidx(intmerit))=0;
                
                if sum(sum(sum(sum(pickr))))==0
                    use_pass_vars=1;
                    cd(design.saveto)
                    load('pass_vars')
                    for roiidxs=1:design.numFolds
                        for meritidxs=1:length(design.merit_thresholds)
                                    for subfold=1:design.numFolds
                                        tmptrackvars(subfold,:)=squeeze(pass_vars(roiidxs, meritidxs, mainfold, subfold,:));
                                    end
                                    varscount=sum(tmptrackvars);
                                    pickr(roiidxs, meritidxs)=not(isempty(find(varscount>=roiidxs)));
                        end
                    end
                end
                
                [roiidx meritidx]=ind2sub([design.numFolds length(design.merit_thresholds)],find(pickr));
                tmpoptimalchk=zeros(design.numFolds,length(roiidx),2);
                for subfold=1:design.numFolds
                    tmpoptimalchk(subfold,find(roiidx==paramidx2use(introi,mainfold,subfold)), introi)=1;
                    tmpoptimalchk(subfold,find(meritidx==paramidx2use(intmerit,mainfold,subfold)), intmerit)=1;
                end
                optimalchk=squeeze(sum(sum(tmpoptimalchk,3)));
                tt=find(optimalchk==max(optimalchk));
                clear diff_params tmpdiff_params
                for diffloop=1:length(tt)
                    for subfold=1:design.numFolds
                        tmpdiff_params(diffloop,subfold)=sum(abs([roiidx(tt(diffloop)) meritidx(tt(diffloop))]-paramidx_perfold(:,subfold)'));
                    end
                    diff_params(diffloop)=sum(tmpdiff_params(diffloop,:));
                end
                paramidxs=tt(find(diff_params==min(diff_params)));
                for subfold=1:design.numFolds
                    paramidx2use(introi,mainfold,subfold)=roiidx(paramidxs(1));
                    paramidx2use(intmerit,mainfold,subfold)=meritidx(paramidxs(1));
                end
                params2use(introi,mainfold)=roiidx(paramidxs(1));
                params2use(intmerit,mainfold)=design.merit_thresholds(meritidx(paramidxs(1)));
                
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
                            trackvars(mainfold,subfold,:)=squeeze(vars2use(mainfold, subfold,tmpidx(introi),tmpidx(intmerit),:));
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
        
        [Beta{mainfold},bint,r,rint,stats] = regress(design.outcome(trainsubjects),[design.data(trainsubjects,Vars2pick_main{mainfold}), design.extradata(trainsubjects,:)]);
        
        GetProbs{mainfold} = glmval([1;Beta{mainfold}],[design.data(testsubjects,Vars2pick_main{mainfold}), design.extradata(testsubjects,:)],design.link);
        design.prediction(testsubjects)=GetProbs{mainfold};
        Merit(mainfold)= -sqrt(abs(design.outcome(testsubjects)-GetProbs{mainfold})'*abs(design.outcome(testsubjects)-GetProbs{mainfold})/length(design.outcome(testsubjects)));
        [r(mainfold), p(mainfold)] = corr(GetProbs{mainfold}, design.outcome(testsubjects));
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
        [prec, tpr, fpr, thresh] = prec_rec_rob_mod_ej(design.prediction', logical(design.outcome),'tstPrecRec', 'plotPR',0);
        fscore=(prec.*tpr)./(prec+tpr);
        stats.F1score=max(fscore);
        save('Results',  'Merit', 'GetProbs', 'Beta', 'params2pick', 'params2use', 'stats', 'Vars2pick_main', 'design');
        disp(['Result for occurrence threshold=' num2str(startROIthresh) ' and merit threshold=' num2str(design.merit_thresholds(startmeritthresh)) ' : AUC=' num2str(stats.overall_AUC)]);
end
end

function [optim_mainfold paramidx_persubfold] = findbestparams_MR(AUCTune_in, design, startROIthresh, startAUCthresh)

for folds=1:design.numFolds*design.numFolds
    [mainfold subfolds]=ind2sub([design.numFolds design.numFolds], folds);
    AUCTune(:,:,mainfold,subfolds)=AUCTune_in(mainfold, subfolds, startROIthresh:end, startAUCthresh:end);
end

ROIthresholds=startROIthresh:design.numFolds;
AUCthresholds=design.merit_thresholds(startAUCthresh:end);
for mainfold=1:10
    for subfold=1:10
        for roi=1:size(AUCTune,1)
            for auc=1:size(AUCTune,2)
                        AUCthreshmatrix(roi,auc,mainfold, subfold)=AUCthresholds(auc);
                        ROIcritmatrix(roi,auc,mainfold, subfold)=ROIthresholds(roi);
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
        critvarperfold=squeeze(AUCTune(:,:,mainfold,subfold));
        AUCthreshperfold=squeeze(AUCthreshmatrix(:,:,mainfold,subfold));
        ROIcritperfold=squeeze(ROIcritmatrix(:,:,mainfold,subfold));
        [bestvalpersubfold(mainfold, subfold) bestparamspersubfold(mainfold,subfold)]=max(critvarperfold(:));
        findmax=find(critvarperfold==max((critvarperfold(:))));
        
            if length(findmax)==1
                winsub(2,subfold)=AUCthreshperfold(findmax);
                winsub(1,subfold)=ROIcritperfold(findmax);
            elseif length(findmax)>1
                best_sub=[];
                best_sub(2,:)=AUCthreshperfold(findmax);
                best_sub(1,:)=ROIcritperfold(findmax);
                
                winsub(:,subfold)=find_optim(best_sub);
                
            else
                winsub(:,subfold)=NaN;
            end

        if length(find(isnan(winsub(1,:))==1))==length(winsub(1,:));
            disp(['Warning: No parameters found for mainfold ' num2str(mainfold) '. Possibly no variables passed the thresholds in this mainfold.'])
        end
    end
    optim_mainfold(:,mainfold)=find_optim(winsub);
    for subfold=1:design.numFolds
        %Stability trheshold index
        paramidx_persubfold(1,mainfold,subfold)=optim_mainfold(1,mainfold);
        %Prediction error threshold index
        clear diffs
        for zz=1:length(design.merit_thresholds)
            diffs(zz)=abs(design.merit_thresholds(zz)-optim_mainfold(2,mainfold));
        end
        ff=find(diffs==min(diffs));
        try
        paramidx_persubfold(2,mainfold,subfold)=ff(1);
        catch ME
            pause
        end
    end
end
end

function winsub=find_optim(bestsub)
findintersectsub=zeros(size(bestsub));
for j=1:2
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
        f=find(isctsub==max(isctsub));
        for maxloop=1:length(f)
            tmp1=bestsub(:,f(maxloop));
            diff(maxloop)=sum(abs(idealval-tmp1'));
        end
        val2use=find(diff==min(diff));
        winsub=bestsub(:,f(val2use(1)));
    end
end