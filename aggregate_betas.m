function [choicefreq sorted_betas betas ALLBetas brainbetas]=aggregate_betas(design, Vars2pick_main, Beta, nbrainFeatures, nCoVars, numFolds, varnames, covarnames);

brainbetas=zeros(nbrainFeatures,numFolds);
pickedbrain=zeros(nbrainFeatures,numFolds);
covarbetas=zeros(nCoVars,numFolds);
pickedcovars=zeros(nCoVars,numFolds);
if isempty(Vars2pick_main)==1
    for mainfold=1:numFolds
        Vars2pick{mainfold}=1:length(varnames);
    end
    Vars2pick_main=Vars2pick;
end
for mainfold=1:numFolds
    for vars=1:length(Vars2pick_main{mainfold})
        brainbetas(Vars2pick_main{mainfold}(vars),mainfold)=Beta{mainfold}(1+vars);
        if Beta{mainfold}(1+vars)~=0
            pickedbrain(Vars2pick_main{mainfold}(vars),mainfold)=1;
        else
            pickedbrain(Vars2pick_main{mainfold}(vars),mainfold)=0;
        end
    end
    if length(Beta{mainfold})>length(varnames)+1
        for covars=1:nCoVars
            covarbetas(covars,mainfold)=Beta{mainfold}(1+length(Vars2pick_main{mainfold})+covars);
            if Beta{mainfold}(1+length(Vars2pick_main{mainfold})+covars)~=0
                pickedcovars(covars,mainfold)=1;
            end
        end
    else
        disp(['Ignored covariates in mainfold ' num2str(mainfold)])
    end
end
choicefreq=sum(pickedbrain,2);
zeromean_brain=mean(brainbetas,2);
brainbetas(brainbetas==0)=NaN;
NaNmean_brain=nanmean(brainbetas,2);
zeromean_cov=mean(covarbetas,2);
covarchoicefreq=sum(pickedcovars,2);
covarbetas(covarbetas==0)=NaN;
NaNmean_cov=nanmean(covarbetas,2);
if length(choicefreq)>length(varnames)
    delthem=1;
    l=length(varnames);
    while delthem==1 && l<length(choicefreq)
        l=l+1;
        if choicefreq(l)~=0
            delthem=0;
        end
    end
    if delthem==1
        ALLBetas=[[choicefreq(1:length(varnames)), zeromean_brain(1:length(varnames)), NaNmean_brain(1:length(varnames))];[covarchoicefreq, zeromean_cov, NaNmean_cov]];
    else
        disp(['something is up with the lengths of varnames and number of actual variables'])
    end
else
    ALLBetas=[[choicefreq, zeromean_brain, NaNmean_brain];[covarchoicefreq, zeromean_cov, NaNmean_cov]];
end
[B idx]=sort(ALLBetas(:,1));
idx=flipud(idx);
ALL_Betas=ALLBetas(idx,:);
BB=[brainbetas;ones(nCoVars,design.numFolds)];
for n=1:size(design.data,2)
    try
    strs{n}=varnames{n};
    catch
        strs{n}=num2str(varnames(n));
    end
end
if length(covarnames)>1
    for n=1:length(covarnames)
        strs{size(design.data,2)+n}=covarnames{n};
    end
elseif isempty(covarnames)==0
    strs{size(design.data,2)+1}=covarnames;
end
sorted_betas(:,1)=strs(idx);
betas(:,1)=strs;
sorted_betas=array2table(sorted_betas);
betas=array2table(betas);
sorted_betas(:,2:4)=array2table(ALL_Betas);
betas(:,2:4)=array2table(ALLBetas);
sorted_betas(:,5:4+design.numFolds)=array2table(BB(idx,:));
tmp={'Variable_index' 'Choice_frequency' 'Mean_beta_incl_0' 'Mean_beta_excl_0'};
for n=1:design.numFolds
    tmp{4+n}=['cv' num2str(n)];
end
sorted_betas.Properties.VariableNames=tmp;
betas.Properties.VariableNames={'Variable_index' 'Choice_frequency' 'Mean_beta_incl_0' 'Mean_beta_excl_0'};
end
