function looper(rootdir, numreps, FS, method);
% this function automatically repeats your analysis a number of times

% example usage:
% looper('/home/lee/local/analysis1', 10, 1, 'EN')

% rootdir: the directory in which your design.mat file is saved and in which the analyses will be saved
% numreps: number of repetitions of the analysis to be run
% FS: whether or not feature selection should be used (1=yes, 0=no)
% method: which regression method should be used (by all_ML), 'EN' (Elatic Net), 'MR' (Multiple Regression), or 'TB' (Treebagger)

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
