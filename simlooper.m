clear

folderstr1='/media/DISK1/Lee/'
for wdh=1
    f_step=[75  200 400  750 1000];
    n_step=[75  200  400  750 1000 2000];
    for n=1:length(n_step)
        for f=1:length(f_step)
            %             disp(['N=' num2str(n_step(n)) '; F=' num2str(f_step(f))]);
            for fs=1:2
                if fs==1
                    fsstr='noFS';
                else
                    fsstr='FS';
                end
                for methodloop=1:3
                    if methodloop==1
                        method2use='EN';
                    elseif methodloop==2
                        method2use='MR';
                    elseif methodloop==3
                        method2use='TB';
                    end
                    for bag=1
                        if bag==1
                            nboot=1;
                            bagstr='nobag';
                        else
                            nboot=25;
                            bagstr='bag';
                        end
                        F=f_step(f);
                        N=n_step(n);
                        for null=1
                            if null==1
                                savestr2use=[folderstr1 filesep fsstr '_' method2use '_' bagstr filesep num2str(F) filesep num2str(N) filesep num2str(wdh)];
                                loadstr=[folderstr1 filesep fsstr '_EN_' bagstr filesep num2str(F) filesep num2str(N) filesep num2str(wdh)];
                            else
                                savestr2use=[folderstr1  filesep fsstr '_' method2use '_' bagstr filesep num2str(F) filesep num2str(N) filesep num2str(wdh)];
                                loadstr=[folderstr1  filesep fsstr '_EN_' bagstr filesep num2str(F) filesep num2str(N) filesep num2str(wdh)];
                            end
                            if exist([savestr2use filesep 'Results.mat'])==0 && exist([savestr2use filesep 'results.mat'])==0
                                if exist([savestr2use filesep 'pass_vars.mat'])==0 && exist([loadstr filesep 'pass_vars.mat'])~=0 && fs==2
                                    load([loadstr filesep 'design.mat']);
                                    load([loadstr filesep 'merit_per_var.mat']);
                                    load([loadstr filesep 'pass_vars.mat']);
                                    design.saveto=savestr2use;
                                    if ~exist(design.saveto, 'dir')
                                        mkdir(design.saveto);
                                    end
                                    save([savestr2use filesep 'design.mat'], 'design');
                                    save([savestr2use filesep 'merit_per_var.mat'], 'merit_per_var');
                                    save([savestr2use filesep 'pass_vars.mat'], 'pass_vars');
                                    stats=all_ML(design, fs-1, method2use);
                                else
                                    if exist([savestr2use filesep 'design.mat'])==0
                                        [design.data design.outcome braincorr outcomecorr]=randomeffects2(1000, 2000, 30, 16, .2, 2.3, [.5 .7 1], [5 5 90], [100], [1], [100], [0.25]);
                                        d1=shuffle(1:size(design.data',1));
                                        d2=shuffle(1:size(design.data',2));
                                        if null==1
                                            design.outcome=shuffle(design.outcome);
                                        end
                                        [design]=create_design(design.data(d2(1:N), d1(1:F))', [1:F], [], [],design.outcome(d2(1:N)), 'linear', nboot, 10,8,8, 'median', 0, savestr2use,[1:N], 'balanced', [], [], 0);
                                    else
                                        load([savestr2use filesep 'design.mat'])
                                        design.saveto=savestr2use;
                                        save([savestr2use filesep 'design.mat'], 'design');
                                    end
                                    stats=all_ML(design, fs-1, method2use);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
