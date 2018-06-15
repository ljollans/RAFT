function [accuracy recall specificity falseposrate falsenegrate precision NPV AUC F1]=class_mets_new(truth, preds)
%

[AUC,fpr,tpr] = fastAUC(logical(truth),preds',0);
[prec, tpr,fpr, thresh] = prec_rec_rob_mod(preds', logical(truth),'tstPrecRec', 'plotPR',0);
fscore=2*((prec.*tpr)./(prec+tpr));
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
nclassaspos=sum(pred);
nclassasneg=length(pred)-nclassaspos;
nposinposclass=sum(truth(find(pred==1)));
nneginnegclass=length(find(pred(find(truth==0))==0));
nposinnegclass=sum(truth(find(pred==0)));
nneginposclass=length(find(pred(find(truth==0))==1));

% of elements in the pos/neg class, what proportion was wrongly classified?
falsenegrate=(nposinnegclass*100)/length(find(truth==1));
falseposrate=(nneginposclass*100)/length(find(truth==0));

% of elements classified as pos/neg, what proportion was actually pos/neg?
precision=(nposinposclass*100)/length(find(pred==1));
NPV=(nneginnegclass*100)/length(find(pred==0));

%of elements in the pos/neg class, what proportion was correctly classed?
recall=(nposinposclass*100)/length(find(truth==1));
specificity=(nneginnegclass*100)/length(find(truth==0));

%overall what proportion of elements was correclty classified?
accuracy=(sum(correct)*100)/length(correct);

end
