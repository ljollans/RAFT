function stats=eliminatereruns(filedir)
% restarts Regularized Adaptive Feature Thresholding where it left off
%
% for comments and questions please contact lee.jollans@gmail.com

savestr=filedir;
cd(filedir)
load('design.mat')
design.saveto=savestr;
if exist([savestr filesep 'Results.mat'])==0
    if exist([savestr filesep 'LHmerit.mat'])==1
        load('LHmerit.mat');
        load('vars2use.mat');
        load('lambda_values.mat');
        load('pass_vars.mat');
        [params2pick, Merit, Beta, Vars2pick_main, GetProbs, design, stats] = RAFT_Model_Selection(design, vars2use, LHmerit, lambda_values, pass_vars, 1, 1, 1);
    elseif exist([savestr filesep 'merit_per_var.mat'])==0
        [params2pick, Merit, Beta, Vars2pick_main, stats]=RAFT(design);
    elseif exist([savestr filesep 'pass_vars.mat'])==0
        load('merit_per_var.mat')
        [design, pass_vars]=RAFT_do_thresh(design, merit_per_var);
        RAFT_2nd_Level(design, pass_vars, [1:(design.numFolds*design.numFolds)]);
        [LHmerit, vars2use, lambda_values] = RAFT_collect_2nd_level(design);
        [params2pick, Merit, Beta, Vars2pick_main, GetProbs, design, stats] = RAFT_Model_Selection(design, vars2use, LHmerit, lambda_values, pass_vars, 1, 1, 1);
    else
        s=dir(filedir);
        finished_folds=[];
        for n=1:size(s,1)
            if length(s(n).name)>20
                if strcmp('Second_level_results', s(n).name(1:20))==1
                    if strcmp('_', s(n).name(23))==1
                        mainfold=str2double(s(n).name(22));
                        if strcmp('.', s(n).name(25))==1
                            subfold=str2double(s(n).name(24));
                        else
                            subfold=str2double(s(n).name(24:25));
                        end
                    else
                        mainfold=str2double(s(n).name(22:23));
                        if strcmp('.', s(n).name(26))==1
                            subfold=str2double(s(n).name(25));
                        else
                            subfold=str2double(s(n).name(25:26));
                        end
                    end
                    tmp=sub2ind([design.numFolds design.numFolds], mainfold, subfold);
                    finished_folds=[finished_folds, tmp];
                end
            end
        end
        folds2do=[];
        for n=1:100
            if isempty(find(finished_folds==n))==1
                folds2do=[folds2do, n];
            end
        end
        load('pass_vars.mat')
        RAFT_2nd_Level(design, pass_vars, folds2do);
        [LHmerit, vars2use, lambda_values] = RAFT_collect_2nd_level(design);
        [params2pick, Merit, Beta, Vars2pick_main, GetProbs, design, stats] = RAFT_Model_Selection(design, vars2use, LHmerit, lambda_values, pass_vars, 1, 1, 1);
    end
else
    load('Results.mat')
end
end