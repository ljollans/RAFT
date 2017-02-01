function looper(rootdir, numreps, FS, method);
for rep=1:numreps;
    for null=1:2
        load([rootdir filesep 'design.mat'])
        if null==1
            design.saveto=[rootdir filesep num2str(rep) '_actual'];
        else
            design.outcome=shuffle(design.outcome);
            design.saveto=[rootdir filesep num2str(rep) '_null'];
        end
        if ~exist(design.saveto, 'dir')
            mkdir(design.saveto);
        end
        if exist([design.saveto filesep 'Results.mat'])==0
            [design.mainfold,design.subfolds]=AssignFolds(size(design.data,1),design.numFolds,design.numFolds);
            chklogfolds=0;
            while chklogfolds==0;
                for n=1:design.numFolds
                    for m=1:design.numFolds
                        c(n,m)=sum(design.outcome(find(design.subfolds(:,n)==m)));
                    end;
                end;
                if isempty(find(c==0))==0
                    [design.mainfold,design.subfolds]=AssignFolds(size(design.data,1),design.numFolds,design.numFolds);
                else
                    chklogfolds=1;
                end
            end
            
            cd(design.saveto);
            save('design', 'design');
            stats=all_ML(design,FS, method);
        end
    end
end
end