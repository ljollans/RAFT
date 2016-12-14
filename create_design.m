function [design]=create_design(data, datalabels, covariates, covarlabels, outcome, type, nboot, numFolds, numLambdas, numAlphas, bagcrit, clean,saveto, subid, balanced, winsorize_data, winsorize_extradata, exclude_binary_vars)
% data: input features
% datalabels: variable names for input features
% covariates: features that will not be included in the feature selection step if this step is performed
% covarlabels: variable names for covariates
% outcome: outcome variable
% type: 'linear' or 'logistic' regression
% nboot: number of bags for bagging, if nboot==1 no bagging is performed
% numFolds: number of cross-validation folds
% numLambdas: number of parameters lambda to evaluate for the elastic net
% numAlphas: number of parameters alpha to evaluate for the elastic net
% bagcrit: criterion used to aggregate bags if bagging is used, 'median' or 'cdf' (cumulative distribution function)
% clean: if this is set to 1 then all datapoints containing NaNs or -100 (missing data) will be removed, as well as features and participants that only have one value
% saveto: location to save results. does not need to be an existing folder, this folder will automatically be created if it does not exist. should be an empty folder
% subid: participant ID codes
% balanced: 'balanced' or 'unbalanced' depending on whether the  groups in logistic regression are balanced or not
% winsorize_data: logical vector indicating for each feature whether it should be winsorized. default is all features are winsorized
% winsorize_extradata: logical vector indicating for each covariate whether it should be winsorized. default is covariates are not winsorized
% exclude_binary_vars: if this is set to 1 then all variables that only have two unique values will be removed
%
% for comments and questions please contact lee.jollans@gmail.com

design.data=data;
if isempty (datalabels)==1
    disp(['please supply variable names'])
end
design.vars=datalabels;
if isempty(covariates)==1
    design.extradata=zeros(length(outcome),1);
    design.covarlabels='X';
else
    design.extradata=covariates;
    if isempty (covarlabels)==1
        disp(['please supply covariate names'])
    end
    design.covarlabels=covarlabels; 
end
design.outcome=outcome; 
design.saveto=saveto;
if isempty (saveto)==1
    disp(['please supply a save location'])
end
if isempty(subid)==1
    design.subid=1:length(outcome);
else
    design.subid=subid;
end
design.nboot=nboot; 
design.numFolds=numFolds; design.nMainFolds=numFolds; design.nSubFolds=numFolds;
design.numLambdas=numLambdas; 
design.numAlphas=numAlphas; 
if strcmp(bagcrit, 'cdf')==1 
    design.bagcrit=bagcrit;
elseif strcmp(bagcrit, 'median')==1
    design.bagcrit=bagcrit;
else
    disp('please select a valid bagging criterion (cdf or mean)')
end

design.siglevel=0.05; %only used for bagcrit='cdf'. 
design.Ratio=2/3; %portion of data used in each bootstrap iteration

design.group=design.outcome;
design.data=design.data';
if length(design.vars)~=size(design.data,2)
        design.data=design.data';
end
design.nvars=size(design.data,2);

if size(design.extradata,1)~=size(design.data,1)
    design.extradata=design.extradata';
end

if clean==1
    s1=sum(design.data,1);     f1=find(s1==0);
    design.data(:,f1)=[];
    if strcmp(class(design.vars), 'cell')==1
        for n=1:length(f1)
            design.vars{f1(n)}=[];
        end
        design.vars=design.vars(~cellfun('isempty', design.vars));
    else
        design.vars(f1)=[];
    end
    design.nvars=size(design.data,2);
    
    s2=sum(design.data,2);f2=find(s2==0);
    design.data(f2,:)=[];design.outcome(f2)=[];design.subid(f2)=[];design.extradata(f2,:)=[];
    
    [r c]=find(isnan(design.data)==1);r=unique(r);
    design.subid(r)=[];design.extradata(r,:)=[];design.data(r,:)=[];design.outcome(r)=[];
    
    [r c]=find(isnan(design.extradata)==1);r=unique(r);
    design.subid(r)=[];design.extradata(r,:)=[];design.data(r,:)=[];design.outcome(r)=[];
    
    [r c]=find(isnan(design.outcome)==1);r=unique(r);
    design.subid(r)=[];design.extradata(r,:)=[];design.data(r,:)=[];design.outcome(r)=[];
    
    [r c]=find((design.data)==-100);r=unique(r);
    design.subid(r)=[];design.extradata(r,:)=[];design.data(r,:)=[];design.outcome(r)=[];
    
    [r c]=find((design.extradata)==-100);r=unique(r);
    design.subid(r)=[];design.extradata(r,:)=[];design.data(r,:)=[];design.outcome(r)=[];
    
    [r c]=find((design.outcome)==-100);r=unique(r);
    design.subid(r)=[];design.extradata(r,:)=[];design.data(r,:)=[];design.outcome(r)=[];
end
if exclude_binary_vars==1
    for n=1:design.nvars
        s3(n)=length(unique(design.data(:,n)));
    end
    f3=find(s3<3);
    design.data(:,f3)=[];
    if strcmp(class(design.vars), 'cell')==1
        for n=1:length(f3)
            design.vars{f3(n)}=[];
        end
        design.vars=design.vars(~cellfun('isempty', design.vars));
    else
        design.vars(f3)=[];
    end
    design.nvars=size(design.data,2);
end
%%
design.group=design.outcome;
if length(design.vars)~=design.nvars
    disp('check length of the variable names vector')
end
if size(design.extradata,1)~=size(design.data,1)
    design.extradata=design.extradata';
end
%%
[design.mainfold,design.subfolds]=AssignFolds(size(design.data,1),design.numFolds,design.numFolds);
design.lambda = logspace(-1,0,design.numLambdas);
design.alpha = linspace(0.01,1.0,design.numAlphas);
if isempty(winsorize_data)==0
    design.winsorize.who.data = winsorize_data; %1 to Winsorize, 0 otherwise
else
    design.winsorize.who.data = ones(1,size(design.data,2)); %1 to Winsorize, 0 otherwise
end
if isempty(winsorize_extradata)==0
    design.winsorize.who.extradata = winsorize_extradata; %1 to Winsorize, 0 otherwise
else
    design.winsorize.who.extradata = zeros(1,size(design.extradata,2)); %1 to Winsorize, 0 otherwise
end
design.winsorize.howmuch=3;
%%
design.type=type;
if strcmp(type, 'linear')==1
    design.merit_thresholds=linspace(-1, -0.8, 10);
    design.distribution='normal';
    design.family='gaussian';
    design.link='identity';
%     design.outcome=zscore(design.outcome);
elseif strcmp(type, 'logistic')==1
    if isempty(balanced)==1
        disp(['please enter "balanced" or "unbalanced"']);;
    else
        design.balanced=balanced;
    end
    design.outcome=logical(design.outcome);
    design.merit_thresholds=linspace(0.51, 0.6, 10);
    design.distribution='binomial';
    design.family='binomial';
    design.link='logit';
    design.posClass2use=1; %make sure design.outcome isn't in logical format
    chklogfolds=0;
    while chklogfolds==0;
        for n=1:design.numFolds
            for m=1:design.numFolds
                c(n,m)=sum(design.outcome(find(design.subfolds(:,n)==m)));
            end;
        end;
        if isempty(find(c==0))==0
            [design.mainfold,design.subfolds]=AssignFolds(size(design.data,1),design.numFolds,design.numFolds);
        else
            chklogfolds=1;
        end
    end
else
    disp('Please enter a valid design type (linear or logistic)');
end
%%
if ~exist(design.saveto, 'dir')
    mkdir(design.saveto);
end

%%

[design.data, mu, sigma] = zscore(design.data);
[r c]=find(design.data>design.winsorize.howmuch);
for n=1:length(r)
    if design.winsorize.who.data(c(n))==1
        design.data(r(n),c(n))  =  design.winsorize.howmuch;
    end
end
[r c]=find(design.data<-design.winsorize.howmuch);
for n=1:length(r)
    if design.winsorize.who.data(c(n))==1
        design.data(r(n),c(n))  =  -design.winsorize.howmuch;
    end
end
[design.extradata, mu, sigma] = zscore(design.extradata);
[r c]=find(design.extradata>design.winsorize.howmuch);
for n=1:length(r)
    if design.winsorize.who.extradata(c(n))==1
        design.extradata(r(n),c(n))  =  design.winsorize.howmuch;
    end
end
[r c]=find(design.extradata<-design.winsorize.howmuch);
for n=1:length(r)
    if design.winsorize.who.extradata(c(n))==1
        design.extradata(r(n),c(n))  =  -design.winsorize.howmuch;
    end
end
design.outcome=zscore(design.outcome);
%%
cd(design.saveto);
save('design', 'design');
end

function [mainfold,subfolds]=AssignFolds(nSubjects, nMainFolds, nSubFolds)
mainfold=crossvalind('Kfold',nSubjects,nMainFolds);
subfolds = -1*ones(nSubjects, nMainFolds);
for subfoldlooper=1:nMainFolds
    SFsubjects = find(mainfold~=subfoldlooper);
    subfolds(SFsubjects,subfoldlooper) = crossvalind('Kfold',length(SFsubjects),nSubFolds);
end
end