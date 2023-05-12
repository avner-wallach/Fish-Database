function [database] = add_experiment()
expstrct=generate_expstructs;
load('z:\analysis_data\ELL_Database.mat');

dlgtitle = ['New Experiment'];
prompt = {'expname','locations'};
dims = [1 35];
definput = {'',''};
answer = inputdlg(prompt,dlgtitle,dims,definput);
locs=str2num(answer{2});
experiment.expname=answer{1};
%check if already exist
if(numel(database))
    K=find(cellfun(@(x) strcmp(x,answer{1}),{database.experiments.expname}));
    if(~numel(K))
        K=numel(database.experiments)+1;
    end
else
    K=1;
end

experiment.locations=repmat(expstrct.locations,numel(locs),1);
for i=1:numel(locs)
    experiment.locations(i).locname=locs(i);
    dlgtitle = ['Experiment ',experiment.expname,' Location ',num2str(locs(i))];
    prompt = {'Groups'};
    dims = [1 100];
    definput = {''};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    grps=str2num(answer{1});
    experiment.locations(i).groups=repmat(expstrct.locations.groups,numel(grps),1);

    for j=1:numel(grps)
        experiment.locations(i).groups(j).groupname=grps(j);
        dlgtitle = ['Experiment ',experiment.expname,' Location ',num2str(locs(i)),' Group ',num2str(grps(j))];
        prompt = {'Wall Seg','Brass Seg'};
        dims = [1 100];
        definput = {'',''};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        experiment.locations(i).groups(j).wallseg=str2num(answer{1});
        experiment.locations(i).groups(j).brassseg=str2num(answer{2});
    end
end
database.experiments(K)=experiment;
save('z:\analysis_data\ELL_Database.mat','database','-append');
end

