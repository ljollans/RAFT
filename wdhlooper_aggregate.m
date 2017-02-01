function [SORTED_ALL_BETAS AUC F1score]=wdhlooper_aggregate(dir2use, numreps, fileprefix, filepostfix, fileappend4save)
% function to aggregate results from several runs of the same model
% example usage: [SORTED_ALL_BETAS AUC F1score]=wdhlooper_aggregate('/home/lee/local/tmp/Results', 5, 'run_', '_null','_5runs_null');

% dir2use is the directory in which the folders for each analysis to be summarized are stored

% numreps is the numebr of models that are to be summarized

% fileprefix and filepostfix are the elements of folder names before and after the integer signalling the number of the model in that folder
% for instance if I name my folders 'run1', 'run2', ... etc then fileprefix is 'run' and filepostfix is empty ('').
% if i name my folders 'repetition_1_nullmodel', 'repetition_2_nullmodel'... etc then fileprefix is 'repetition_' and filepostfix is '_nullmodel'

% fileappend4save is the string appended to .csv files with results that are saved in dir2use. I use this to keep track of what analysis 
% which .csv file belongs to and how many models are summarized

% SORTED_ALL_BETAS: table with table with variable names, in how many CV folds they were chosen on average, 
% mean beta value if folds where this variable was not chosen are taken as beta=0, 
% and mean beta value if folds in which this variable was not chosen were ignored
% variables are sorted by beta value. This Table is also saved as a .csv file automatically.

% AUC: vector containing all AUC (if logistic regression) or R (if lienar regression) values for all summarized analyses

% F1score: vector containing all F1 scores (if logistic regression) or P values (if lienar regression) values for all summarized analyses

thresholdsused=[];
thresholdspicked=[];
for wdh=1:numreps
    clear stats r p design
    load([dir2use filesep fileprefix num2str(wdh) filepostfix filesep 'Results.mat']);
    load([dir2use filesep fileprefix num2str(wdh) filepostfix filesep 'design.mat']);
    switch design.type
        case 'logistic'
            if exist('stats', 'var')==1
                AUC(wdh)=stats.overall_AUC;
                F1score(wdh)=stats.F1score;
            else
                %whatever thi is without feature selection
            end
        case 'linear'
            if exist('stats', 'var')==1
                AUC(wdh)=stats.r;
                F1score(wdh)=stats.p;
            else
                AUC(wdh)=r;
                F1score(wdh)=p;
            end
    end
    if exist('Vars2pick_main', 'var')==0
        Vars2pick_main=[];
    end
    [choicefreq{wdh} sorted_Betas{wdh} Betas{wdh} Unsorted_allbetas{wdh}]=aggregate_betas(Vars2pick_main, Beta, size(design.data,2), size(design.extradata,2), design.numFolds, design.vars, design.covarlabels);
    choicefreq_all(:,wdh)=Unsorted_allbetas{wdh}(:,1);
    Meanbeta_incl0(:,wdh)=Unsorted_allbetas{wdh}(:,2);
    Meanbeta_excl0(:,wdh)=Unsorted_allbetas{wdh}(:,3);
    if exist('params2use', 'var')==1
        thresholdsused=[thresholdsused, params2use];
        thresholdspicked=[thresholdspicked, params2pick];
    else
        thresholdsused=[thresholdsused; [alpha2use', lambda2use']];
    end
end
if exist('params2use', 'var')==1
    thresholds=array2table([thresholdsused; thresholdspicked]');
    thresholds.Properties.VariableNames={'occurrence_threshold_used', 'merit_thresholds_used', 'lambda_used', 'alpha_used', 'occurrence_threshold_picked', 'merit_thresholds_picked', 'lambda_picked', 'alpha_picked'};
else
    thresholds=array2table([thresholdsused]);
    thresholds.Properties.VariableNames={'alpha_used', 'lambda_used'};
end
writetable(thresholds, [dir2use filesep 'Thresholds' fileappend4save '.xls']);
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
X=array2table([AUC', F1score']);
switch design.type
    case 'logistic'
        X.Properties.VariableNames={'AUC', 'F1score'};
    case 'linear'
        X.Properties.VariableNames={'R', 'P'};
end
writetable(X, [dir2use filesep 'Modelfit' fileappend4save '.xls']) 
end
