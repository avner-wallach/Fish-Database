function [ strct ] = get_regression_models( varargin )
%GET_POLAR_MODEL Calculate polar decoding model parameters

pxl2mm=str2num(getenv('PXLSIZE'));
%% params
% params.mode='cart'; %cart=cartesian ; polar=polar
params.hidmode='BP'; %BP=backprop; GC=granule-cell (random); only in NI models
% params.nhidden=100;
% params.ngc=1e3; %number of 'granule cells'
params.freqmem=0;
params.freqmode='on';
params.t_col=[1 5:14];     %input columns
% params.ref='circ'; %circ- circular tank; rect- rectangular tank; pole- pole object
params.rmax=inf; %max distance from object to use for modeling
params.sampnum=5e4;
params.priorbins=50;
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
% params.lfb_cols=z_col;
params.fb_cols=find(cellfun(@(x) numel(strfind(x,'lfp')),data_struct(1).fnames));

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
x=[]; y=[]; a=[]; z=[]; t=[]; fb=[];% lfb=[];
for k=1:params.K
    x1=data_struct(k).data(:,params.xy_cols(1));
    y1=data_struct(k).data(:,params.xy_cols(2));
    a1=data_struct(k).data(:,params.a_col);
    
    x=[x;x1];
    y=[y;y1];
    a=[a;a1];
    z=[z;data_struct(k).data(:,z_col)];
    if(params.t_col~=0)
        t=[t;data_struct(k).data(:,params.t_col)];
        if(params.freqmem>0)
            iei1=t(:,iiei);
            for i=1:params.freqmem
                iei1=[nan;iei1(1:end-1)];
                t=[iei1 t];
                fnames={['iei',num2str(i)],fnames{:}};
            end
            iiei=find(cellfun(@(x) numel(strfind(x,'iei')),fnames));
        end
        
    end
    if(numel(params.fb_cols) & params.fb_cols~=0)
        fb=[fb;data_struct(k).data(:,params.fb_cols)];
    end
            
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

%% location grid
M=params.priorbins;
p=1;
loc=[x y];
x_edges=linspace(min(x)*p,max(x)*p,M);
y_edges=linspace(min(y)*p,max(y)*p,M+1);
x_bins=edge2bin(x_edges);
y_bins=edge2bin(y_edges);
loc_bins={x_bins,y_bins};
%% past input
% z=medfilt1(z,3);
% fb=medfilt1(fb,3);
fbp=[nan(1,size(fb,2));fb(1:end-1,:)]; %past sensory input
tp=[nan(1,size(t,2));t(1:end-1,:)]; %past motor input
[ net_predictor,rsq_predictor,sig_predictor] = ni_network_new([t tp fbp]',z','BP');
z_pred=net_predictor([t tp fbp]')';
%% downsample
if(~isinf(params.sampnum) & numel(z)>params.sampnum)
    ind=randperm(numel(z),params.sampnum);
    z=z(ind);
    loc=loc(ind,:);
    t=t(ind,:);
    fb=fb(ind,:);
    z_pred=z_pred(ind);
%     fbp=fbp(ind,:);
end
%% models

% ex-afferent regrssion 
[ model(1).net,model(1).rsq,model(1).sig] = ni_network_new(loc',z','BP');
z_ex=model(1).net(loc'); %ex-afference signal
model(1).R=nancorr(z_ex(:),z(:))^2;
model(1).R_ex=model(1).R;
model(1).name='EX';
model(1).cols=[];

z_re=z'-z_ex;   %reafference component

% names={'CD','LFB','FB','CD*LFB','CD+FB','CD*FB','DC*FBp'};
names={'CD','LFB','FB','CD*LFB','CD*FB','CDLOC','CDLFBp'};
inputs={t,z,fb,[t z],[t fb],[t loc],[t z_pred]};
dims={[],[],[],[],[],[],[]};
cols={params.t_col,z_col,params.fb_cols,[params.t_col z_col],[params.t_col params.fb_cols],[nan],[nan]};
for i=1:numel(names)
    model(i+1)=get_model(inputs{i},z_re,'GC',dims{i},names{i},cols{i});
end
%% single electrode addition
[fbcols,I]=setdiff(params.fb_cols,z_col);
[g2lmodel.net,g2lmodel.rsq,g2lmodel.sig]=ni_network_new(fb(:,I)',z','BP');
for i=1:numel(fbcols)    
    [l2gmodel(i).net,l2gmodel(i).rsq,l2gmodel(i).sig]=ni_network_new(z',fb(:,I(i))','BP');
%     smodel(i)=get_model([z fb(:,I(i))],z_re,'GC',[],['LFB+FB',num2str(i)],[z_col fbcols(i)]);
%     scdmodel(i)=get_model([t z fb(:,I(i))],z_re,'GC',[],['CD+LFB+FB',num2str(i)],[params.t_col z_col fbcols(i)]);
    c(i)=nancorr(z,fb(:,I(i)));
end

%%
%prior pdf p(x,y)
[X,Y]=meshgrid(loc_bins{1},loc_bins{2});
[f,xi] = ksdensity(loc,[X(:) Y(:)]);
F=log(f);
F(isinf(F))=min(F(~isinf(F)));
[ net_prior,rsq_prior,sig_prior] = ni_network_new([X(:) Y(:)]',F','BP');

%normalization pdf p(z)
xi=linspace(nanmin(z),nanmax(z),1e4);
[f,xi] = ksdensity(z,xi(:));
[ net_norm,rsq_norm,sig_norm ] =ni_network_new(xi(:)',f(:)','BP',1,25);

%% output
strct.z_col=z_col;
strct.fb_col=params.fb_cols;
strct.t_col=params.t_col;
strct.model=model;
strct.l2gmodel=l2gmodel;
strct.g2lmodel=g2lmodel;
strct.c=c;

strct.norm.net=net_norm;
strct.norm.rsq=rsq_norm;
strct.norm.sig=sig_norm;

strct.prior.net=net_prior;
strct.prior.rsq=rsq_prior;
strct.prior.sig=sig_prior;

strct.loc.max=max(loc,[],1);
strct.loc.min=min(loc,[],1);
strct.loc.mean=[nanmean(loc(:,1)) nanmean(loc(:,2))];

strct.mot.max=max(t,[],1);
strct.mot.min=min(t,[],1);
strct.mot.mean=nanmean(t,1);

strct.decoding=[];
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

    function model=get_model(inputs,targets,mode,dims,name,cols)
        [ model.net,model.rsq,model.sig ] = ni_network_new(inputs',targets,mode,dims);
        [ ex_net,ex_rsq,ex_sig] = ni_network_new(loc',(z'-model.net(inputs')),'BP');
        model.R=nancorr(z_ex(:),z-model.net(inputs')')^2;
        model.R_ex=ex_rsq;        
        model.name=name;
        model.cols=cols;
    end

end

