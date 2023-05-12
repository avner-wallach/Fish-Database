function [ decoding_struct ] = evaluate_regression_models( varargin )
%GET_POLAR_MODEL Calculate polar decoding model parameters

pxl2mm=str2num(getenv('PXLSIZE'));
figure;
COL=colormap('lines');
% COL(4,:)=0.5*[1 1 1];
fcolor=[0 0 0];
% fcolor=[1 1 1];
close(gcf);
%% params
params.mode='cart'; %cart=cartesian ; polar=polar
params.nhidden=100;
params.cornersize=180;
params.freqmem=0;
params.freqmode='on';
params.fontsize=14;
params.msize=8;
params.t_col=[1 5:14];     %input columns
params.loc_nbins=25;
params.minsamp=100;
params.plotting=0;
params.bgcol=[1 1 1];
%% get varins
n=numel(varargin);
Narg=3;
if(n<Narg | n>Narg&mod(n-Narg,2))
    error('too few input arguments!');
end
data_struct=varargin{1};
file_struct=varargin{2};
model_struct=varargin{3};

fb_col=model_struct.fb_col;
t_col=model_struct.t_col;
z_col=model_struct.z_col;

params.xy_cols=[find(cellfun(@(x) (strcmp(x,'X')),data_struct(1).fnames)) ...
    find(cellfun(@(x) (strcmp(x,'Y')),data_struct(1).fnames))];%default is spatial mapping
params.a_col=find(cellfun(@(x) (strcmp(x,'azim')),data_struct(1).fnames));

if(n>3)
    i=4;
    while(i<n)
        params.(varargin{i})=varargin{i+1};
        i=i+2;
    end    
end
params.K=numel(data_struct);
fnames=data_struct(1).fnames(params.t_col);
iiei=find(cellfun(@(x) numel(strfind(x,'iei')),fnames));
%% get data vectors
x=[]; y=[]; a=[]; z=[]; d=[];%t=[]; fb=[];% lfb=[];
for k=1:params.K
    x1=data_struct(k).data(:,params.xy_cols(1));
    y1=data_struct(k).data(:,params.xy_cols(2));
    a1=data_struct(k).data(:,params.a_col);
    
    x=[x;x1];
    y=[y;y1];
    a=[a;a1];
    z=[z;data_struct(k).data(:,z_col)];
    d=[d;data_struct(k).data];
            
end
x=x*pxl2mm;
y=y*pxl2mm;
   
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
loc=[x y];
%% location grid
M=params.loc_nbins;
p=1;
x_edges=linspace(min(x)*p,max(x)*p,M);
y_edges=linspace(min(y)*p,max(y)*p,M+1);
loc_edges={x_edges,y_edges};
x_bins=edge2bin(x_edges);
y_bins=edge2bin(y_edges);
loc_bins={x_bins,y_bins};

%% evaluate mutual info using regression model
% prior=exp(model_struct.prior.net(loc')');
% normz=model_struct.norm.net(z')';
[x,ex]=discretize(x,x_edges);
[y,ey]=discretize(y,y_edges);
[tmp1,tmp2,xy]=unique(1000*x+y); %convert to single discrete vector
xy(isnan(x+y))=nan;
Hxy=entropy(xy);

ind=randperm(numel(z));
zz=z(ind);
[N0,Uxy(:,:,1),U(1),Pxy(:,:,1)]=get_info(); %control
names{1}='CTL';

zz=z;
[N0,Uxy(:,:,2),U(2),Pxy(:,:,2)]=get_info(); %no NI
names{2}='EX';

for k=2:numel(model_struct.model)
    inputs=d(:,model_struct.model(k).cols)';
    z_ni=model_struct.model(k).net(inputs)';
    zz=z-z_ni;
    [N0,Uxy(:,:,k+1),U(k+1),Pxy(:,:,k+1)]=get_info();
    names{k+1}=model_struct.model(k).name;
end

if(params.plotting)
    Fall=figure;
    K=k+1;
    r1=ceil(sqrt(K));
    [ha, pos] = tight_subplot(r1, r1, [.03 .01],[.02 .02],[.02 .02]);
    for i=1:K        
            axes(ha(i));
            plot_surface(Uxy(:,:,i),N0);
            colormap(flipud(colormap('pink')));
            set(gca,'CLim',[0 1]);
            title([names{i},'; U=',num2str(U(i))]);
    end
end

decoding_struct.names=names;
decoding_struct.N=N0;
decoding_struct.Uxy=Uxy;
decoding_struct.Pxy=Pxy;
decoding_struct.Hxy=Hxy;
decoding_struct.U=U;
decoding_struct.Umean=permute(nanmean(Uxy,[1 2]),[2 3 1]);
decoding_struct.Umedian=permute(nanmedian(Uxy,[1 2]),[2 3 1]);
return;

%% evaluate decoding using regression model
prior=exp(model_struct.prior.net(loc')');
normz=model_struct.norm.net(z')';

z_ex=model_struct.model(1).net(loc')';
ind=randperm(numel(z));
zz=z(ind);
normzz=normz(ind);
[N0,S(:,:,1),L(:,:,1),R(:,:,1),P(:,:,1)]=get_posterior(); %control
sig(1)=nanstd(zz(:)-z_ex(:));
zz=z;
normzz=normz;
[N0,S(:,:,2),L(:,:,2),R(:,:,2),P(:,:,2)]=get_posterior(); %no NI
sig(2)=nanstd(zz(:)-z_ex(:));


for k=2:numel(model_struct.model)
    if(~strcmp(model_struct.model(k).name,'DC*FBp'))
        inputs=d(:,model_struct.model(k).cols)';
    else
        inputs=cdfbp';
    end
    z_ni=model_struct.model(k).net(inputs)';
    zz=z-z_ni;
    zi=linspace(nanmin(zz),nanmax(zz),1e3);
    [f,zi] = ksdensity(zz,zi(:));
    normzz=interp1(zi,f,zz);  
    
    sig(k+1)=nanstd(zz(:)-z_ex(:));
    [N0,S(:,:,k+1),L(:,:,k+1),R(:,:,k+1),P(:,:,k+1)]=get_posterior();

end
names={'CTL' model_struct.model.name};
if(params.plotting)
    Fall=figure;
    [ha, pos] = tight_subplot(4, k+1, [.03 .01],[.02 .02],[.02 .02]);
    MM={S,L,R,1-P}; 
    rows={'S','L','R','P'};
    clim={[0 2],[0 3],[1 2.5],[.95 1]};
    k=1;
    for i=1:numel(MM)        
        for j=1:size(MM{i},3)
            axes(ha(k));
            plot_surface(MM{i}(:,:,j),N0);
            colormap(flipud(colormap('pink')));
            set(gca,'CLim',clim{i});
            if(i==1)
                title(names{j});
            end
            if(j==1)
                ylabel(rows{i});
            end                
            k=k+1;
        end
    end
end

%% output
decoding_struct.names=names;
decoding_struct.N=N0;
decoding_struct.S=S;
decoding_struct.L=L;
decoding_struct.R=R;
decoding_struct.P=P;
decoding_struct.sig=sig;

%% functions

    function plot_surface(Mz,A)
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
        surf(x_bins,y_bins,zeros(size(Mz')),...
            'CData',Mz',...
            'LineStyle','none','FaceAlpha','interp','FaceColor','interp'...
            ,'AlphaData',A','AlphaDataMapping','scaled');
        % a=max(quantile(N0(:),params.maxa),10);
        set(gca,'ALim',[50 75]);
        view(0,90);        
        set(gca,'YDir','normal');
        axis('image');
        set(gca,'Xlim',x_bins([1 end]),'Ylim',y_bins([1 end]));
        hold off;
        cb=colorbar;
        set(gca,'Color','none','XColor','none','YColor','none');
        set(gcf,'Color',params.bgcol);
        set(cb,'Color',1-params.bgcol,'Location','east','AxisLocation','out');
%         cb.Position = cb.Position + [0.015 0 -.01 0];
%         cb.Ticks = [-2 -1 0 1 2];

    end

    function [N0,Uxy,U,Pxy]=get_info()
        Mx=numel(loc_edges{1})-1;
        My=numel(loc_edges{2})-1;
        Bn=12;
        Bq=[0.01 0.99];
        N0=zeros(Mx,My);
        Uxy=nan(Mx,My);
        Pxy=nan(Mx,My);
        U=nan;
        [pz,zedges]=histcounts(zz,Bn,'Normalization','probability','BinLimits',quantile(zz,Bq));
        pz=pz/sum(pz); %normalize
        Hz=nansum(-pz.*log2(pz));
%         [Hz,Sz] = entropy(zz);
        
        for i=1:Mx
            for j=1:My
                idx=find(inrange(loc(:,1),[loc_edges{1}(i) loc_edges{1}(i+1)]) &...
                         inrange(loc(:,2),[loc_edges{2}(j) loc_edges{2}(j+1)]));
                if(numel(idx)>params.minsamp)                                        
                    pzz=histcounts(zz(idx),zedges,'Normalization','probability');
                    pzz=pzz/sum(pzz); %normalize
                    Hzz=nansum(-pzz.*log2(pzz));
%                     PZZ{i,j}=pzz; Vzz(i,j)=nanvar(zz(idx)); Mzz(i,j)=nanmean(zz(idx));
                    Uxy(i,j)=(Hz-Hzz);
                    Pxy(i,j)=nansum(pzz.^2./pz);
                    N0(i,j)=numel(idx);                
                end
            end
        end
        U=nansum(Uxy.*N0,'all')/sum(N0(:));        
    end
        
    function [N0,S,L,R,P]=get_posterior()
        Mx=numel(loc_edges{1})-1;
        My=numel(loc_edges{2})-1;
        N=1e3;
        N0=zeros(Mx,My);
        S=nan(Mx,My);
        L=S;
        R=S;        
        pval=nan(Mx,My);
%         pmean=zeros(Mx,My);

        for i=1:Mx
            for j=1:My
                idx=find(inrange(loc(:,1),[loc_edges{1}(i) loc_edges{1}(i+1)]) &...
                         inrange(loc(:,2),[loc_edges{2}(j) loc_edges{2}(j+1)]));
                if(numel(idx)>params.minsamp)                                        
                    M_EX(i,j)=nanmean(zz(idx));
%                     M_EX(i,j)=nanmean(z_ex(idx));
                    S(i,j)=nanstd(zz(idx)-M_EX(i,j));
                    l=1/sqrt(2*pi*S(i,j)^2)*exp(-(zz(idx)-M_EX(i,j)).^2/(2*S(i,j)^2));
                    r=l./normzz(idx);                   
                    L(i,j)=nanmean(l);
                    R(i,j)=nanmean(r);

                    if(sum(~isnan(r)))
%                         pval(i,j)=ranksum(r,r_ctl,'tail','right');
                        pval(i,j)=signrank(r,1,'tail','right');
                    end
                end
                N0(i,j)=numel(idx);
            end
        end        
        P=pval;
%         P=(R.*N0)./sum(N0(:));
    end

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


end

