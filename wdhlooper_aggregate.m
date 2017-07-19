function [ALL_BETAS SORTED_ALL_BETAS AUC F1score]=wdhlooper_aggregate(dir2use, numreps, fileprefix, filepostfix, fileappend4save)

%clear
%dir='/media/hanni/707E536D7E532ADC/Google Drive/MS ML/SPM/parcellated_modeling/P_out_P3bVismo0/Results_CMSMLlinCOGNz_P3bVis';
%numreps=8;
thresholdsused=[];
thresholdspicked=[];
for wdh=1:numreps
    clear stats r p design
    load([dir2use filesep fileprefix num2str(wdh) filepostfix filesep 'design.mat']);
    load([dir2use filesep fileprefix num2str(wdh) filepostfix filesep 'Results.mat']);
    if exist([dir2use filesep fileprefix num2str(wdh) filepostfix filesep 'prediction.mat'])
        design.prediction=pred;
    end

    switch design.type
        case 'logistic'
            if exist('stats', 'var')==1
                x=stats.tpr-stats.fpr;
                pred=zeros(size(design.prediction)); 
                try
                    t=find(x==max(x));
                    pred(find(design.prediction>=stats.thresh(t(1))))=1;
                catch
                    disp('oops...')
                end
                for n=1:length(pred), correct(n)=isequal(pred(n), design.outcome(n)); end
                accuracy(wdh)=(sum(correct)*100)/length(correct);
                sensitivity(wdh)=(sum(correct(find(design.outcome==1)))*100)/length(correct(find(design.outcome==1)));
                specificity(wdh)=(sum(correct(find(design.outcome==0)))*100)/length(correct(find(design.outcome==0)));
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
    if isempty(design.covarlabels)==1
        for n=1:size(design.extradata,2)
        design.covarlabels{n}=num2str(n);
        end
    end
    [choicefreq{wdh} sorted_Betas{wdh} Betas{wdh} Unsorted_allbetas{wdh} brainbetas{wdh}]=aggregate_betas(design, Vars2pick_main, Beta, size(design.data,2), size(design.extradata,2), design.numFolds, design.vars, design.covarlabels);
    choicefreq_all(:,wdh)=Unsorted_allbetas{wdh}(:,1);
    Meanbeta_incl0(:,wdh)=Unsorted_allbetas{wdh}(:,2);
    Meanbeta_excl0(:,wdh)=Unsorted_allbetas{wdh}(:,3);
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
        X=array2table([AUC', F1score', accuracy', sensitivity', specificity']);
        X.Properties.VariableNames={'AUC', 'F1score', 'accuracy', 'sensitivity', 'specificity'};
    case 'linear'
        X=array2table([AUC', F1score']);
        X.Properties.VariableNames={'R', 'P'};
end
writetable(X, [dir2use filesep 'Modelfit' fileappend4save '.xls']) 
end
