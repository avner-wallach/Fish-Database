function [ group ] = process_group(eod0,file0,ops,group)
%PROCESS_GROUP get all stats for group

tcols=[6 1];
T=numel(tcols);
setenv('PXLSIZE','0.423'); %pxl to mm conversion

ttitles={'Tail','Rate'};
groupnum=group.groupname;
fnames={'units','lfp','lfpratio','wall','obj'};
for i=1:numel(fnames)
    if(isfield(group,fnames{i}))
        group=rmfield(group,fnames{i});
    end
end

if(numel(group.brassseg))
    segnum=group.brassseg;
    objidx=1;
else
    segnum=group.wallseg;
    objidx=0;
end
% mode='wall';
mode='';
process;

% mode='obj';
% process;

    function process
    if(~numel(segnum))
        group.lfp=[];
        group.units=[];
        return;
    end

    eod=eod0(segnum);
    file=file0(segnum);

    %get cols
    lfpcol=find(cellfun(@(x) numel(strfind(x,['lfp',num2str(groupnum)])),eod.fnames));
    allunitcols=find(cellfun(@(x) numel(strfind(x,['sc'])),eod.fnames));
    indlim_inds=find(cellfun(@(x) numel(strfind(x,['sc',num2str(groupnum)])),eod.fnames(allunitcols)));
    unitcols=allunitcols(indlim_inds);
    U=numel(unitcols);
    % unitcols=find(cellfun(@(x) numel(strfind(x,['sc',num2str(groupnum)])),eod.fnames));
    ind_lim=ops.seg(segnum).ind_lim(indlim_inds,:);


    %LFP measures
    tr=nanmedian(eod.avtraces(:,:,ops.seg(segnum).LFPgroups_pca{groupnum},:),4);
    [m,mind]=min(tr,[],2);
    [M,Mind]=max(tr,[],2);
    r=(M-(abs(m)))./(M+abs(m));
    group.lfp.([mode,'ratio'])=mean(median(r,1));


    % response map entropy
    [Mz,N0]=plot_tuning2d(eod,file,lfpcol,'mfunc','mean','objidx',objidx,'plotting',0);
    [group.lfp.([mode,'Hx']),group.lfp.([mode,'Sx'])]=entropy(Mz);

    for j=1:T 
        tcol=tcols(j);
        [Mz,N0]=plot_tuning2d(eod,file,lfpcol,'t_col',tcol,'mfunc','slope','objidx',objidx,'plotting',0);
        group.lfp.([mode,'C_',ttitles{j}])=nansum(abs(Mz).*N0,'all')/nansum(N0(:));

        [ Izt,Izt_xy,Izxy,Izxy_t,izt_xy ] = compute_info(eod,file,lfpcol,'t_col',tcol,'objidx',objidx); 
        group.lfp.([mode,'U_',ttitles{j}])=Izt; %uncrtainty coef. for variable
        group.lfp.([mode,'Uxy_',ttitles{j}])=Izt_xy; %uncrtainty coef. for variable, given position
    end

    %% process group units
    minisi=0.05; %minimal isi before EOD to compute baseline fr
    baseline_bin=[-0.02 -0.005];        
    max_edges=[0:1e-3:0.04];   

%     [db,ustrct]=generate_expstructs();
    indu=find(ops.seg(segnum).Spkgroups==groupnum);

%     group=rmfield(group,'units');
    for i=1:numel(unitcols)
        if(isfield(group,'units'))
            if(numel(group.units)>=i)
                ustrct=group.units(i);
            end
        end
        ustrct.uname=eod.fnames{unitcols(i)};
%         %use xcor to find spiking mode
%         bins=file.units.xcor.bins;
%         val=file.units.xcor.val{indu}(:,i,i);
%         [M,Mind]=max(val(bins>0));
%         bpos=bins(bins>0);
%         ustrct.([mode,'spmode'])=bpos(Mind);
        ustrct.([mode,'spmode'])=file.units.isi.mode{indu}(i);
        ustrct.([mode,'refractory'])=file.units.isi.refractory{indu}(i);
        ustrct.([mode,'burst'])=file.units.isi.burst{indu}(i);

        %use raster to compute baseline and peak response
        ind=find(eod.data(:,1)>minisi);
        ind=ind(inrange(ind,ind_lim(i,:)));
        spind=find(ismember(eod.raster{i}(:,2),ind));
        ustrct.([mode,'baseline'])=sum(inrange(eod.raster{i}(spind,1),baseline_bin))/diff(baseline_bin)/numel(ind);
        h=histcounts(eod.raster{i}(spind,1),max_edges)/mean(diff(max_edges))/numel(ind);
        ustrct.([mode,'max'])=max(h);

        % response map entropy
        [Mz,N0]=plot_tuning2d(eod,file,unitcols(i),'mfunc','mean','objidx',objidx,...
            'ind_lim',ind_lim(i,:),'plotting',0);
        [ustrct.([mode,'Hx']),ustrct.([mode,'Sx'])]=entropy(Mz);

        ttcols=[tcols lfpcol];
        tttitles=[ttitles 'Lfp'];
        for j=1:numel(ttcols) 
            tcol=ttcols(j);
            [Mz,N0]=plot_tuning2d(eod,file,unitcols(i),'t_col',tcol,'mfunc','slope','objidx',objidx,...
                'ind_lim',ind_lim(i,:),'plotting',0);
            ustrct.([mode,'C_',tttitles{j}])=nansum(abs(Mz).*N0,'all')/nansum(N0(:));

            [ Izt,Izt_xy,Izxy,Izxy_t,izt_xy ] = compute_info(eod,file,unitcols(i),'t_col',tcol,...
                'ind_lim',ind_lim(i,:),'objidx',objidx); 
            ustrct.([mode,'U_',tttitles{j}])=Izt; %uncrtainty coef. for variable
            ustrct.([mode,'Uxy_',tttitles{j}])=Izt_xy; %uncrtainty coef. for variable, given position
        end
        Us(i)=ustrct;
    end
    group.units=Us;    
    end
end
