function [b, pred, r, p, mse]=do_mreg(X, truth, nboot, saveto)

kb=[ones(size(X,1),1), X];
[mf, sf]=AssignFolds(size(X,1),10,10);
%%mreg
clear btmp
if nboot>1
    for boot=1:nboot
        parfor n=1:10
            ftrain=find(mf~=n);
            [Xboot,Yboot,indexselect]=bootstrapal(kb(ftrain,:),truth(ftrain),2/3);
            [btmp(boot,:,n),bint,r,rint,stats] = regress(Yboot,Xboot);
        end
    end
    b=squeeze(mean(btmp,1));
else
    parfor n=1:10
        ftrain=find(mf~=n);
        [b(:,n),bint,r,rint,stats] = regress(truth(ftrain), kb(ftrain,:));
    end
end
for n=1:10
    ftest=find(mf==n);
    pred(ftest)=glmval([b(:,n)],kb(ftest,2:end), 'identity');
end
[r p]=corr(pred', truth);
mse=(abs(truth-pred')'*abs(truth-pred')/length(truth));

if ~exist(saveto, 'dir')
    mkdir(saveto);
end
stats.r=r;
stats.p=p;
stats.mse=mse;
cd(saveto)
save('data', 'X', 'truth');
save('betas', 'b');
save('prediction', 'pred');
save('Results',  'pred', 'b', 'r', 'p', 'mse');
end