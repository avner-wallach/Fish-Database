%% process all experiments
% load('z:\analysis_data\ELL__Models_Database.mat');
mode=0; %1=brass, 0=wall
i=1;
tic
while(i<=numel(database.experiments))    
    database.experiments(i).expname    
    database=process_experiment_models(database,database.experiments(i).expname,mode);
    save('z:\analysis_data\ELL_Models_Database_wall.mat','database','-append');
    toc
    i=i+1;
end
% return;
%% get recording length
for i=1:numel(database.experiments)
    load(['z:\analysis_data\',num2str(database.experiments(i).expname),'\output.mat'],'eod');
    recdur(i)=sum(arrayfun(@(x) nansum(x.data(:,1)),eod)/3600);
end
return;
%% collect data
groups=[];
for i=1:numel(database.experiments)
%     r=[];
    for k=1:numel(database.experiments(i).groups)        
        G=database.experiments(i).groups(k);
        gname=[database.experiments(i).expname,'_',num2str(k)];
%         h=1/G.decoding.Hxy;
%         T=[cell2table({gname},'VariableNames',{'gname'}),cell2table({G.model.R},'VariableNames',strcat({G.model.name},'_Model')),...            
%             array2table(G.decoding.U*h,'VariableNames',strcat(G.decoding.names,'_Decode')),...
%             cell2table({G.model.rsq},'VariableNames',strcat({G.model.name},'_RSQ')),...            
%             array2table(mean(abs(G.c)),'VariableNames',{'C'}),...
%             array2table(numel(G.c),'VariableNames',{'nFB'}),...
%             array2table(G.g2lmodel.rsq,'VariableNames',{'G2L'}),...
%             array2table(mean(cell2mat({G.l2gmodel.rsq})),'VariableNames',{'L2G'}),...
%             array2table(G.tail_v,'VariableNames',{'Tail_Var'}),...
%             array2table(G.tail_v_ctl,'VariableNames',{'Tail_Var_Ctl'}),...
%             array2table(G.tail_Uxy,'VariableNames',{'Tail_Uxy'}),...
%             array2table(G.tail_Uxy_ctl,'VariableNames',{'Tail_Uxy_Ctl'}),...
%             array2table(G.rate_v,'VariableNames',{'Rate_Var'}),...
%             array2table(G.rate_v_ctl,'VariableNames',{'Rate_Var_Ctl'}),...
%             array2table(G.rate_Uxy,'VariableNames',{'Rate_Uxy'}),...
%             array2table(G.rate_Uxy_ctl,'VariableNames',{'Rate_Uxy_Ctl'}),...
%             array2table(G.C_total,'VariableNames',strcat({G.model([1 2 3 5 6]).name},'_Reafference'))];

        v_none=nanvar(G.M_re(:).*G.N0(:)/nansum(G.N0(:)));
        for p=1:numel(G.model)
            v(p)=nanvar(G.model(p).M_re(:).*G.N0(:)/nansum(G.N0(:)));
        end
        T=[cell2table({gname},'VariableNames',{'gname'}),...
            cell2table({G.model.rsq},'VariableNames',strcat({G.model.name},'_RSQ')),...            
            array2table(mean(abs(G.c)),'VariableNames',{'Cmean'}),...
            array2table(max(abs(G.c(G.c<.999))),'VariableNames',{'Cmax'}),...
            array2table(numel(G.c),'VariableNames',{'nFB'}),...
            cell2table({G.R_ex G.model.R_ex},'VariableNames',strcat({'NONE' G.model.name},'_R_EX')),...                        
            array2table(G.tail_v_flat,'VariableNames',{'Var_Flat'}),...
            array2table(G.tail_v_flat_ctl,'VariableNames',{'Var_Flat_Ctl'}),...
            array2table(G.tail_v_weighted,'VariableNames',{'Var_Weight'}),...
            array2table(G.tail_v_weighted_ctl,'VariableNames',{'Var_Weight_Ctl'}),...
            array2table([v_none v],'VariableNames',strcat({'NONE' G.model.name},'_Var')),...
            array2table(G.tail_Uxy,'VariableNames',{'Tail_Uxy'}),...
            array2table(G.tail_Uxy_ctl,'VariableNames',{'Tail_Uxy_Ctl'}),...
            array2table([G.C_flat cell2mat({G.model.C_flat})],'VariableNames',strcat({'NONE' G.model.name},'_ReFlat')),...
            array2table([G.C_weighted cell2mat({G.model.C_weighted})],'VariableNames',strcat({'NONE' G.model.name},'_ReWeight'))];
        
        groups=[groups;T];
    end
%     r0=repmat(r(:,5)-r(:,1),1,k);   
%     ind=ones(size(r0))-eye(size(r0));    
%     R0=[R0;r0(find(ind))];
%     R_all=[R_all;r];
end
save('ELL__Models_Database.mat','groups','-append');
%% display stats to select examples
T=[groups.FB_RSQ groups.FBp_RSQ groups.CD_RSQ groups.CDFB_RSQ groups.CDFBp_RSQ groups.CDLOC_RSQ];
T2=T(:,1:end-1); p=size(T2,2);
T2=T2./(T(:,end)*ones(1,p));

P=[groups.NONE_ReFlat groups.FB_ReFlat groups.FBp_ReFlat groups.CD_ReFlat groups.CDFB_ReFlat groups.CDFBp_ReFlat groups.CDLOC_ReFlat];
P2=P(:,2:end); p=size(P2,2);
P2=P2./(P(:,1)*ones(1,p));

Q=[groups.NONE_ReWeight groups.FB_ReWeight groups.FBp_ReWeight groups.CD_ReWeight groups.CDFB_ReWeight groups.CDFBp_ReWeight groups.CDLOC_ReWeight];
Q2=Q(:,2:end); p=size(Q2,2);
q=Q(:,1)-Q(:,end);
Q2=Q2./(Q(:,1)*ones(1,p));
% Q2=(Q2-Q(:,end))./(q*ones(1,p));

V=[groups.NONE_Var groups.FB_Var groups.FBp_Var groups.CD_Var groups.CDFB_Var groups.CDFBp_Var groups.CDLOC_Var];
V2=V(:,2:end); p=size(V2,2);
V2=V2./(V(:,1)*ones(1,p));
% 
% 
% T=[groups.EX_Model+2e-3*randn(size(groups.EX_Model)) groups.CD_Model groups.LFB_Model groups.("CD*LFB_Model") groups.("CD*FB_Model")];
% TT=[T nan(size(T,1),1)];
% TTT=reshape(TT',[],1);
% X=1:numel(TT);
% P=[groups.EX_Decode groups.CD_Decode groups.LFB_Decode groups.("CD*LFB_Decode") groups.("CD*FB_Decode")];
PP=[P2(:,[3 5 6]) nan(size(P2,1),1)];
X=1:numel(PP);
PPP=reshape(PP',[],1);
figure;
% yyaxis('left')
% plot(X,TTT,':*');
% hold on;
% Vline(X(isnan(TTT)));
% yyaxis('right')
plot(X,PPP,':*');
Vline(X(isnan(PPP)));
Xt=X(isnan(PPP))-3;
set(gca,'XTick',Xt,'XTickLabel',{groups.gname{:}},'XTickLabelRotation',90,'XLim',[0 max(X)]);


QQ=[Q2(:,[3 5 6]) nan(size(Q2,1),1)];
X=1:numel(QQ);
QQQ=reshape(QQ',[],1);
figure;
plot(X,QQQ,':*');
Vline(X(isnan(QQQ)));
Xt=X(isnan(QQQ))-3;
set(gca,'XTick',Xt,'XTickLabel',{groups.gname{:}},'XTickLabelRotation',90,'XLim',[0 max(X)]);


