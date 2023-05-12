% all experiments
dbfile='z:\analysis_data\ELL_Database.mat';
load(dbfile);
for i=1:numel(database.experiments)
    expname=database.experiments(i).expname
    database=process_experiment(database,expname);
    save(dbfile,'database','-append');
end

% collect all data
lfp=[];
units=[];
for i=1:numel(database.experiments)
    exp=database.experiments(i);
    exp.expname
    for j=1:numel(exp.locations)
        loc=exp.locations(j);
        for k=1:numel(loc.groups)
            gr=loc.groups(k);
            grname={[exp.expname,'_',num2str(loc.locname),'_',num2str(gr.groupname)]};
            if(numel(gr.lfp))
                lfp=[lfp;[table(grname,'VariableNames',{'grname'}),struct2table(gr.lfp)]];
            end
            if(numel(gr.units))
                for u=1:numel(gr.units)
                    units=[units;[cell2table({grname{1},gr.lfp.ratio},'VariableNames',{'grname','ratio'})...
                        ,struct2table(gr.units(u))]];
                end
            end            
        end
    end
end