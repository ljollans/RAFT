function stats=all_ML(design, feature_selection, method)
% selects correct code to run Elastic Net/ Multiple Regression
% with/without feature selection
%
% for comments and questions please contact lee.jollans@gmail.com

% latest update: may 30th 2017

disp(['Doing ' method ' with FS=' num2str(feature_selection)]);
cd(design.saveto);
if feature_selection==1 & strcmp(method, 'EN')==1
    stats=eliminatereruns(design.saveto)
elseif feature_selection==1 & strcmp(method, 'MR')==1
stats=eliminatererunsMR(design.saveto)
elseif feature_selection==0 & strcmp(method, 'EN')==1
    [Beta, pred, stats]=do_EN(design.data, design.outcome, design.nboot, design.bagcrit, design.saveto, design.type);
elseif feature_selection==0 & strcmp(method, 'MR')==1
    [stats.b, stats.pred, stats.r, stats.p, stats.mse]=do_mreg(design.data, design.outcome, design.nboot, design.saveto);
elseif feature_selection==0 & strcmp(method, 'TB')==1
    [stats]=do_random_forest_matlab(design);
elseif feature_selection==1 & strcmp(method, 'TB')==1
    if exist([design.saveto filesep 'merit_per_var.mat'])==0
        if exist([design.saveto filesep 'tmpmerit.mat'])
            load([design.saveto filesep 'tmpmerit.mat'])
            [merit_per_var] = RAFT_FS(design, tmpmerit);
        else
            [merit_per_var] = RAFT_FS(design, []);
        end
    else
        load([design.saveto filesep 'merit_per_var.mat']);
    end
    if exist([design.saveto filesep 'pass_vars.mat'])==0
        [design, pass_vars]=RAFT_do_thresh(design, merit_per_var);
    else
        load([design.saveto filesep 'pass_vars.mat']);
    end
    if exist([design.saveto filesep 'Merit.mat'])==0
        folds2do=getfolds2do(design.saveto);
         Treebagger_2nd_Level(design, pass_vars, folds2do);
         [Merit]=Treebagger_collect_2nd_level(design);
    else
        load([design.saveto filesep 'Merit.mat']);
    end
    [params2pick, Merit, Vars2pick_main, GetProbs, design, stats] = Treebagger_Model_Selection(design, Merit, pass_vars);
else
    disp('did not do analysis')
    stats=0;
end
end
