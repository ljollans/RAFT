function [choicefreq sorted_betas betas ALLBetas]=aggregate_betas(Vars2pick_main, Beta, nbrainFeatures, nCoVars, numFolds, varnames, covarnames);
% This function gives you some summary statistics for your model.
% load 'Results.mat' into the workspace and run, example:

%[choicefreq sorted_betas betas ALLBetas]=aggregate_betas(Vars2pick_main, Beta, length(design.vars), ...
% length(design.covarlabels), design.numFolds, design.vars, design.covarlabels);

%choicefreq: vector showing in how many CV folds each variable was chosen

% betas: table with variables, in how many CV folds they were chosen, 
% mean beta value if folds where this variable was not chosen are taken as beta=0, 
% and mean beta value if folds in which this variable was not chosen were ignored

% sorted_betas: the same as betas but sorted by mean beta value

% ALLBetas: matrix with the last three columns of sorted_betas

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
        pickedbrain(Vars2pick_main{mainfold}(vars),mainfold)=1;
    end
    %if length(Beta{mainfold})==length(varnames)+length(covarnames)+1
    for covars=1:nCoVars
        covarbetas(covars,mainfold)=Beta{mainfold}(1+length(Vars2pick_main{mainfold})+covars);
        if Beta{mainfold}(1+length(Vars2pick_main{mainfold})+covars)~=0
            pickedcovars(covars,mainfold)=1;
        end
    end
    %else
    %    disp(['Ignored covariates in mainfold ' num2str(mainfold)])
    %end
end
choicefreq=sum(pickedbrain,2);
zeromean_brain=mean(brainbetas,2);
brainbetas(brainbetas==0)=NaN;
NaNmean_brain=nanmean(brainbetas,2);
zeromean_cov=mean(covarbetas,2);
covarchoicefreq=sum(pickedcovars,2);
covarbetas(covarbetas==0)=NaN;
NaNmean_cov=nanmean(covarbetas,2);
ALLBetas=[[choicefreq, zeromean_brain, NaNmean_brain];[covarchoicefreq, zeromean_cov, NaNmean_cov]];
[B idx]=sort(ALLBetas(:,1));
idx=flipud(idx);
ALL_Betas=ALLBetas(idx,:);
for n=1:length(varnames)
    strs{n}=varnames{n};
end
if length(covarnames)>1
    for n=1:length(covarnames)
        strs{length(varnames)+n}=covarnames{n};
    end
else
    strs{length(varnames)+1}=covarnames;
end
sorted_betas(:,1)=strs(idx);
betas(:,1)=strs;
sorted_betas=array2table(sorted_betas);
betas=array2table(betas);
sorted_betas(:,2:4)=array2table(ALL_Betas);
betas(:,2:4)=array2table(ALLBetas);
sorted_betas.Properties.VariableNames={'Variable_index' 'Choice_frequency' 'Mean_beta_incl_0' 'Mean_beta_excl_0'};
betas.Properties.VariableNames={'Variable_index' 'Choice_frequency' 'Mean_beta_incl_0' 'Mean_beta_excl_0'};
end
