function [ALL_BETAS, SORTED_ALL_BETAS, stats1]=wdhlooper_aggregate(dir2use, numreps, fileprefix, filepostfix, fileappend4save)

%clear
%dir='/media/hanni/707E536D7E532ADC/Google Drive/MS ML/SPM/parcellated_modeling/P_out_P3bVismo0/Results_CMSMLlinCOGNz_P3bVis';
%numreps=8;
thresholdsused=[];
thresholdspicked=[];
whatfields={'vars','covarlabels','extradata','numFolds'};
for wdh=1:numreps
    clear stats r p design pred
    for assignfields=1:length(whatfields)
        eval(['design.' whatfields{assignfields} '=[];']);
    end
    design.data=[];
    design.outcome=[];
    load([dir2use filesep fileprefix num2str(wdh) filepostfix filesep 'design.mat']);
    %% artefact from MS ML data where Lee saved design files as xcomp for some reason (why did i do that? i don't know)
    if exist('xcomp') && isempty(design.data)
        design=xcomp;
    end
    %%
    try
        load([dir2use filesep fileprefix num2str(wdh) filepostfix filesep 'Results.mat']);
    catch
        load([dir2use filesep fileprefix num2str(wdh) filepostfix filesep 'results.mat']);
    end
    if exist([dir2use filesep fileprefix num2str(wdh) filepostfix filesep 'prediction.mat'])
        load([dir2use filesep fileprefix num2str(wdh) filepostfix filesep 'prediction.mat']);
        design.prediction=pred;
    end
    try
        load([dir2use filesep fileprefix num2str(wdh) filepostfix filesep 'data.mat']);
        if isempty(design.data),
            design.data=X;
        end;
        if isempty(design.outcome),
            design.outcome=truth;
        end
        assignauto={[1:size(design.data,2)],{1},{NaN(size(design.outcome))},{10}};
        for assignfields=1:length(whatfields)
            if isempty(eval(['design.' whatfields{assignfields}]))
                eval(['design.' whatfields{assignfields} '=' assignauto{assignfields} ';']);
            end
        end
        if isempty(design.type)
            if length(unique(truth))<3
                design.type='logistic';
            else
                design.type='linear';
            end
        end
    end
    switch design.type
        case 'logistic'
                [stats1.accuracy(wdh) stats1.trueposrate(wdh) stats1.truenegrate(wdh) stats1.falseposrate(wdh) stats1.falsenegrate(wdh) stats1.posretrieval(wdh) stats1.negretrieval(wdh) stats1.AUC(wdh) stats1.F1score(wdh)]=class_mets_new(design.outcome, design.prediction);
            
        case 'linear'
            if exist('stats', 'var')==1
                stats1.r(wdh)=stats.r;
                stats1.p(wdh)=stats.p;
            else
                stats1.r(wdh)=r;
                stats1.p(wdh)=p;
            end
    end
    if exist('Vars2pick_main', 'var')==0
        Vars2pick_main=[];
    end
    if isempty(design.covarlabels)==1
        for n=1:size(design.extradata,2)
            design.covarlabels{n}=num2str(n);
        end
    end
    [choicefreq{wdh} sorted_Betas{wdh} Betas{wdh} Unsorted_allbetas{wdh} brainbetas{wdh}]=aggregate_betas(design, Vars2pick_main, Beta, size(design.data,2), size(design.extradata,2), design.numFolds, design.vars, design.covarlabels);
try
    choicefreq_all(:,wdh)=Unsorted_allbetas{wdh}(:,1);
    Meanbeta_incl0(:,wdh)=Unsorted_allbetas{wdh}(:,2);
    Meanbeta_excl0(:,wdh)=Unsorted_allbetas{wdh}(:,3);
catch
    pause
end
end

mcf=mean(choicefreq_all');
mbi=mean(Meanbeta_incl0');
mbe=nanmean(Meanbeta_excl0');
ALL_BETAS(:,1)=Betas{1}(:,1);
ALL_BETAS(:,2)=array2table(mcf');
ALL_BETAS(:,3)=array2table(mbi');
ALL_BETAS(:,4)=array2table(mbe');
ALL_BETAS.Properties.VariableNames={'Variable_index' 'MEANChoice_frequency' 'Mean_beta_incl_0' 'Mean_beta_excl_0'};
[B idx]=sort(table2array(ALL_BETAS(:,2)));
idx=flipud(idx);
SORTED_ALL_BETAS(:,1)=Betas{1}(idx,1);
SORTED_ALL_BETAS(:,2)=array2table(mcf(idx)');
SORTED_ALL_BETAS(:,3)=array2table(mbi(idx)');
SORTED_ALL_BETAS(:,4)=array2table(mbe(idx)');
SORTED_ALL_BETAS.Properties.VariableNames={'Variable_index' 'MEANChoice_frequency' 'Mean_beta_incl_0' 'Mean_beta_excl_0'};
%
writetable(SORTED_ALL_BETAS, [dir2use filesep 'Sorted_Beta_Values' fileappend4save '.xls'])
writetable(ALL_BETAS, [dir2use filesep 'Unsorted_Beta_Values' fileappend4save '.xls'])
%
switch design.type
    case 'logistic'
        X=array2table([stats1.AUC', stats1.F1score', stats1.accuracy', stats1.trueposrate', stats1.truenegrate', stats1.falseposrate', stats1.falsenegrate']);
        X.Properties.VariableNames={'AUC', 'F1score', 'accuracy', 'truepos', 'trueneg', 'falsepos','falseneg'};
    case 'linear'
        X=array2table([stats1.r', stats1.p']);
        X.Properties.VariableNames={'R', 'P'};
end
writetable(X, [dir2use filesep 'Modelfit' fileappend4save '.xls'])
save([dir2use filesep 'Modelfit' fileappend4save '.m'], 'X');
end
