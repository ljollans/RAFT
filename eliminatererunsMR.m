function stats=eliminatererunsMR(filedir)
% restarts Regularized Adaptive Feature Thresholding where it left off
%
% for comments and questions please contact lee.jollans@gmail.com

% latest update: may 30th 2017

savestr=filedir;
cd(filedir)
load('design.mat')
design.saveto=savestr;
if exist([savestr filesep 'Results.mat'])==0
disp('1')
    if exist([savestr filesep 'LHmerit.mat'], 'file') | exist([savestr filesep 'LHmerit10.mat'], 'file')
disp('2')
         [params2pick, Merit, Beta, Vars2pick_main, GetProbs, design, stats] = MR_3rd_level(design,1);
    elseif exist([savestr filesep 'merit_per_var.mat'])==0
        if exist([design.saveto filesep 'tmpmerit.mat'], 'file')
disp('4')
            load([design.saveto filesep 'tmpmerit.mat'])
            [merit_per_var] = RAFT_FS(design, tmpmerit);
            [design, pass_vars]=RAFT_do_thresh(design, merit_per_var);
            MR_2nd_level(design, pass_vars, [1:(design.numFolds*design.numFolds)]);
            if size(design.data,2)>=500
                MR_collect_2nd_level_memorysave(design);
            else
                MR_collect_2nd_level(design);
            end
           
        else
disp('5')
            [merit_per_var] = RAFT_FS(design, []);
            [design, pass_vars]=RAFT_do_thresh(design, merit_per_var);
            MR_2nd_level(design, pass_vars, [1:(design.numFolds*design.numFolds)]);
            if size(design.data,2)>=500
                MR_collect_2nd_level_memorysave(design);
            else
                MR_collect_2nd_level(design);
            end
            [params2pick, Merit, Beta, Vars2pick_main, GetProbs, design, stats] = MR_3rd_level(design, 1, 1);
        end
    elseif exist([savestr filesep 'pass_vars.mat'])==0
disp('6')
        load('merit_per_var.mat')
        [design, pass_vars]=RAFT_do_thresh(design, merit_per_var);
                    MR_2nd_level(design, pass_vars, [1:(design.numFolds*design.numFolds)]);
            if size(design.data,2)>=500
                MR_collect_2nd_level_memorysave(design);
            else
                MR_collect_2nd_level(design);
            end
            [params2pick, Merit, Beta, Vars2pick_main, GetProbs, design, stats] = MR_3rd_level(design, 1, 1);
    else
disp('7')
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
load('design.mat')
        load('pass_vars.mat')

            MR_2nd_level(design, pass_vars, folds2do);

            if size(design.data,2)>=500
                MR_collect_2nd_level_memorysave(design);
            else
                MR_collect_2nd_level(design);
            end
            [params2pick, Merit, Beta, Vars2pick_main, GetProbs, design, stats] = MR_3rd_level(design,1);
    end
else
disp('8')
    load('Results.mat')
end
end
