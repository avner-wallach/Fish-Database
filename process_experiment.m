function [database] = process_experiment(database,expname)
K=find(cellfun(@(x) strcmp(x,expname),{database.experiments.expname}));
if(~numel(K))
    return;
end

%load data
load(['z:\analysis_data\',expname,filesep,'output.mat']);

for i=1:numel(database.experiments(K).locations)
    for j=1:numel(database.experiments(K).locations(i).groups)
        groups(j)=process_group(eod,file,ops,database.experiments(K).locations(i).groups(j));
    end
    database.experiments(K).locations(i).groups=groups;
end

end

