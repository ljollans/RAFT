function [accuracy sensitivity specificity AUC F1]=class_mets(truth, preds)

[AUC,fpr,tpr] = fastAUC(logical(truth),preds',0);
[prec, tpr,fpr, thresh] = prec_rec_rob_mod(preds', logical(truth),'tstPrecRec', 'plotPR',0);
fscore=(prec.*tpr)./(prec+tpr);
F1=max(fscore);

x=tpr-fpr;
pred=zeros(size(preds));
try
    t=find(x==max(x));
    pred(find(preds>=thresh(t(1))))=1;
catch
    disp('oops...')
end
for n=1:length(pred),
    correct(n)=isequal(pred(n), truth(n));
end
accuracy=(sum(correct)*100)/length(correct);
sensitivity=(sum(correct(find(truth==1)))*100)/length(correct(find(truth==1)));
specificity=(sum(correct(find(truth==0)))*100)/length(correct(find(truth==0)));
end