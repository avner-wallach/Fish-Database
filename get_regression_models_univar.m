function [ strct ] = get_regression_models_univar( varargin )
%GET_POLAR_MODEL Calculate polar decoding model parameters

pxl2mm=0.423;%str2num(getenv('PXLSIZE'));
%% params
params.t_col=[6];     %input column
params.sampnum=2e5;
params.minsamp=10;
params.x_nbins=25;
params.y_nbins=25;
params.bgcol=[1 1 1];
params.upsamp=1;
params.plotting=0;
params.nettype='BP';

%% get varins
n=numel(varargin);
if(n<3 | n>3&~mod(n,2))
    error('too few input arguments!');
end
data_struct=varargin{1};
file_struct=varargin{2};
z_col=varargin{3};  %output variable

params.xy_cols=[find(cellfun(@(x) (strcmp(x,'X')),data_struct(1).fnames)) ...
    find(cellfun(@(x) (strcmp(x,'Y')),data_struct(1).fnames))];%default is spatial mapping
params.a_col=find(cellfun(@(x) (strcmp(x,'azim')),data_struct(1).fnames));
% params.fb_cols=setdiff(find(cellfun(@(x) numel(strfind(x,'lfp')),data_struct(1).fnames)),z_col);
params.fb_cols=find(cellfun(@(x) numel(strfind(x,'lfp')),data_struct(1).fnames));

if(n>3)
    i=4;
    while(i<n)
        params.(varargin{i})=varargin{i+1};
        i=i+2;
    end    
end
params.K=numel(data_struct);
%% get data vectors
x=[]; y=[]; a=[]; z=[]; t=[]; fb=[];
for k=1:params.K
    x1=data_struct(k).data(:,params.xy_cols(1));
    y1=data_struct(k).data(:,params.xy_cols(2));
    a1=data_struct(k).data(:,params.a_col);
    
    x=[x;x1];
    y=[y;y1];
    a=[a;a1];
    z=[z;data_struct(k).data(:,z_col)];
    t=[t;data_struct(k).data(:,params.t_col)];        
    
    if(numel(params.fb_cols) & params.fb_cols~=0)
        fb=[fb;data_struct(k).data(:,params.fb_cols)];
    end
            
end
x=x*pxl2mm;
y=y*pxl2mm;
   
%standartize
z=(z-nanmean(z))/nanstd(z);
t=(t-nanmean(t))/nanstd(t);

%% default object idx
if(~isfield(params,'objidx'))
    params.objidx=numel(file_struct(1).objects);
end
%% 
if(params.objidx==0)
    convert_to_wall;
else
    convert_to_object;
end

%% location grid
p=1;
% loc=[x y];
params.x_edges=linspace(min(x)*p,max(x)*p,params.x_nbins+1);
params.y_edges=linspace(min(y)*p,max(y)*p,params.y_nbins+1);
params.x_bins=edge2bin(params.x_edges);
params.y_bins=edge2bin(params.y_edges);
% loc_bins={params.x_bins,params.y_bins};
%% past input
fbp=[nan(1,size(fb,2));fb(1:end-1,:)]; %past sensory input
tp=[nan(1,size(t,2));t(1:end-1,:)]; %past motor input
% [ net_predictor,rsq_predictor,sig_predictor] = ni_network_new([t tp fbp]',fb','BP');
% fb_pred=net_predictor([t tp fbp]')';
%% downsample
if(~isinf(params.sampnum) & numel(z)>params.sampnum)
    ind=randperm(numel(z),params.sampnum);
    z=z(ind);
    x=x(ind);
    y=y(ind);
    t=t(ind,:);
    fb=fb(ind,:);
%     z_pred=z_pred(ind);
    fbp=fbp(ind,:);
end

%% get ex, re maps, vectors
zz=z;
[N0,M_re,M_ex,Pv1]=fitting2d();
z_ex=interp2(params.x_bins,params.y_bins,M_ex',x,y)';
z_re=interp2(params.x_bins,params.y_bins,M_re',x,y)';
C_flat=nanmean(abs(M_re(:)));
C_weighted=nansum(abs(M_re).*N0,'all')/nansum(N0(~isnan(M_re)));
R_ex=nancorr(z_ex(:),z(:))^2;

if(params.plotting)
    F=figure;
    set(gcf,'Units','normalized','Position',[0 0 1 1]);
    [ha,pos]=tight_subplot(1,7,[0 0],[0 0 ],[0 0],[],[]);
    axes(ha(1));
    plot_map(M_re,Pv1);
end

%% models
names={'CD','CDLOC','FB','CDFB','FBp','CDFBp'};
inputs={t,[t x y],fb,[t fb],fbp,[t fbp]};
target=z_re.*t';
% target=z'-z_ex;

for k=1:numel(names)
    model(k)=get_model(inputs{k},target,params.nettype,names{k});
    if(params.plotting)
        axes(ha(k+1));
        plot_map(model(k).M_re,model(k).Pv);
    end
end

for i=1:numel(params.fb_cols)    
    c(i)=nancorr(z,fb(:,i));
end

%% output
strct.z_col=z_col;
strct.fb_col=params.fb_cols;
strct.t_col=params.t_col;
strct.model=model;
strct.c=c;
strct.C_flat=C_flat;
strct.C_weighted=C_weighted;
strct.R_ex=R_ex;
strct.M_ex=M_ex;
strct.M_re=M_re;
strct.Pv=Pv1;
strct.N0=N0;
strct.x_bins=params.x_bins;
strct.y_bins=params.y_bins;

% strct.loc.max=max(loc,[],1);
% strct.loc.min=min(loc,[],1);
% strct.loc.mean=[nanmean(loc(:,1)) nanmean(loc(:,2))];
% 
% strct.mot.max=max(t,[],1);
% strct.mot.min=min(t,[],1);
% strct.mot.mean=nanmean(t,1);

% strct.decoding=[];
%% coordinate transforms
    function convert_to_wall()
        file_struct(1).circle = file_struct(1).circle*pxl2mm;
        rx=file_struct(1).circle(3)/2; ry=file_struct(1).circle(4)/2;    %ellipse radii
        x0=file_struct(1).circle(1) + rx; y0=file_struct(1).circle(2) + ry; %center point
        phi=atan2((y-y0),(x-x0));    %azimuth in tank
        R1=hypot((x-x0),(y-y0)); %distance of fish from center
        R2=(rx*ry)./sqrt((ry*cos(phi)).^2 + (rx*sin(phi)).^2); %distance of nearest point from cetner
        ra=R2-R1;    %distance of fish to nearest wall
        ra(ra<0)=nan;
        th=phi-a;   %egocentric angle of closest wall
        x=ra.*sin(th);
        y=ra.*cos(th);
    end

    function convert_to_object()
          objx=file_struct(k).objects(params.objidx).x*pxl2mm;
          objy=file_struct(k).objects(params.objidx).y*pxl2mm;
          [x,y]=allo2ego(objx,objy,a,x,y); %obj coordinates rel. to LED
    end

    function model=get_model(inputs,targets,mode,name)
        [ model.net,model.rsq,model.sig ] = ni_network_new(inputs',targets,mode,[]);
        zz=z-model.net(inputs')';
        [n0,m_re,m_ex,pv1]=fitting2d();
        model.M_ex=m_ex;
        model.M_re=m_re;
        model.Pv=pv1;

%         [ ex_net,ex_rsq,ex_sig] = ni_network_new(loc',(z'-model.net(inputs')),'BP');
        model.R_ex=nancorr(z_ex(:),z-model.net(inputs')')^2;
        model.C_flat=nanmean(abs(model.M_re(:)));
        model.C_weighted=nansum(abs(model.M_re).*n0,'all')/nansum(n0(~isnan(model.M_re)));
        model.name=name;
    end

    function [N0,P1,P2,Pv1]=fitting2d()
        N0=zeros(params.x_nbins,params.y_nbins);
        P1=nan(params.x_nbins,params.y_nbins);
        P2=nan(params.x_nbins,params.y_nbins);
        Pv1=nan(params.x_nbins,params.y_nbins);
        for i=1:params.x_nbins
            for j=1:params.y_nbins
                idx=find(x>params.x_edges(i) & x<=params.x_edges(i+1) &...
                         y>params.y_edges(j) & y<=params.y_edges(j+1));
                idx=idx(~isnan(t(idx)+zz(idx)));             
                N0(i,j)=numel(idx);
                if(numel(idx)>=max(params.minsamp,3))
                    [f,g]=fit(t(idx),zz(idx),'poly1');
                    P1(i,j)=f.p1;
                    P2(i,j)=f.p2;
                    [rp,pval1]=corr(t(idx),zz(idx));
                    Pv1(i,j)=-log10(pval1);
                end
            end
        end
    end

    function plot_map(Mz,AL)
        if(params.upsamp>1)
            [Y,X]=meshgrid(params.y_bins,params.x_bins);
            F=scatteredInterpolant(X(~isnan(Mz)),Y(~isnan(Mz)),Mz(~isnan(Mz)),'linear','none');
            x_bins=linspace(params.x_bins(1),params.x_bins(end),params.x_nbins*params.upsamp);
            y_bins=linspace(params.y_bins(1),params.y_bins(end),params.y_nbins*params.upsamp);
            [Yq,Xq]=meshgrid(y_bins,x_bins);
            Mz=F(Xq,Yq);    
            N0(isnan(N0))=0;
            AL(isnan(AL))=0;
            N0=interp2(Y,X,N0,Yq,Xq);
            AL=interp2(Y,X,AL,Yq,Xq);
        else
            x_bins=params.x_bins;
            y_bins=params.y_bins;
        end
        
        if(isfield(params,'image'))
            sX=size(params.image,2);
            sY=size(params.image,1);
            ix=([1:sX]-sX/2)*pxl2mm;
            iy=(-[1:sY]+sY/3)*pxl2mm;
            if(params.bgcol(1)==0)
                imagesc(ix,iy,255-params.image);  
            else
                imagesc(ix,iy,params.image);  
            end
        end

        hold on;
        S=surf(x_bins*pxl2mm,y_bins*pxl2mm,zeros(size(Mz')),...
            'CData',Mz',...
            'LineStyle','none','FaceAlpha','interp','FaceColor','interp'...
            ,'AlphaData',AL','AlphaDataMapping','scaled');
        view(0,90);
        set(gca,'YDir','normal');
        axis('image');
        set(gca,'Xlim',x_bins([1 end])*pxl2mm,'Ylim',y_bins([1 end])*pxl2mm);
        set(gca,'Color','none','XColor','none','YColor','none');
        set(gcf,'Color',params.bgcol);
        if(params.bgcol==[1 1 1])
            colormap(gca,flipud(brewermap(64,'RdBU')));
        else
            colormap(gca,invert_map(flipud(brewermap(64,'RdBU'))));
        end
        set(gca,'alim',[1 2.1]);    
        if(isfield(params,'clim'))
            set(gca,'CLim',params.clim);
        else
            set(gca,'Clim',[-1 1]);
        end
    end


end

