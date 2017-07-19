function [mf sf]=check_CV(mf, sf,outcome)
okactual=0;
while okactual==0
    who=0;
    ok=1;
    [mf, sf]=AssignFolds(length(outcome),length(unique(mf)), size(sf,2));
    while ok==1;
        if who<length(unique(mf))*size(sf,2)
            who=who+1;
            [mainfold subfolds]=ind2sub([length(unique(mf)) size(sf,2)], ok);
            train = find(mf ~= mainfold & sf(:,mainfold) ~=subfolds);
            test = find(mf ~= mainfold & sf(:,mainfold)==subfolds);
            if length(unique(outcome(train)))==1
                ok=0;
                okactual=1;
            end
        else
            ok=0;
            okactual=1;
        end
    end
end
            