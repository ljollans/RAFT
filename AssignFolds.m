function [mainfold,subfolds]=AssignFolds(nSubjects, nMainFolds, nSubFolds)
mainfold=crossvalind('Kfold',nSubjects,nMainFolds);
subfolds = -1*ones(nSubjects, nMainFolds);
for subfoldlooper=1:nMainFolds
    SFsubjects = find(mainfold~=subfoldlooper);
    subfolds(SFsubjects,subfoldlooper) = crossvalind('Kfold',length(SFsubjects),nSubFolds);
end
end