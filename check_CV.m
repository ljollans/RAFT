function [mf sf]=check_CV(mf, sf,outcome, numFolds)
okactual=0;
while okactual==0
    who=0;
    ok=1;
    [mf, sf]=AssignFolds(size(outcome,1),numFolds,numFolds);
    while ok==1;
        if who<length(unique(mf))*size(sf,2)
            who=who+1;
            [mainfold subfolds]=ind2sub([length(unique(mf)) size(sf,2)], who);
            clear train test
            train = find(mf ~= mainfold & sf(:,mainfold) ~=subfolds); ltr(who)=length(unique(outcome(train)));
            test = find(mf ~= mainfold & sf(:,mainfold)==subfolds); lte(who)=length(unique(outcome(test)));
            if length(unique(outcome(train)))==1 | length(unique(outcome(test)))==1
                ok=0;
                okactual=0;
            end
        else
            ok=0;
            okactual=1;
        end
    end
end
            