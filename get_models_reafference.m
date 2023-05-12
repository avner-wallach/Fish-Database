function [ C_total ] = get_models_reafference( varargin )
%GET_POLAR_MODEL Calculate polar decoding model parameters

%% params
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

if(n>3)
    i=4;
    while(i<n)
        params.(varargin{i})=varargin{i+1};
        i=i+2;
    end    
end
%% get data vectors
z=data_struct.data(:,z_col);
cd=data_struct.data(:,t_col);
fb=data_struct.data(:,fb_col);
z_cd=nanzscore(z-model_struct.model(2).net(cd')');
z_lfb=nanzscore(z-model_struct.model(3).net(z')');
z_cdlfb=nanzscore(z-model_struct.model(5).net([cd z]')');
z_cdfb=nanzscore(z-model_struct.model(6).net([cd fb]')');

data_struct.data(:,z_col)=nanzscore(z);
data_struct.data(:,z_col+1)=z_cd;   data_struct.fnames{z_col+1}='CD';
data_struct.data(:,z_col+2)=z_lfb;  data_struct.fnames{z_col+2}='LFB';
data_struct.data(:,z_col+3)=z_cdlfb;data_struct.fnames{z_col+3}='CDLFB';
data_struct.data(:,z_col+4)=z_cdfb; data_struct.fnames{z_col+4}='CDFB';

for i=1:5
    [Mz,N0]=plot_tuning2d(data_struct,file_struct,z_col+i-1,'mfunc','slope','t_col',z_col,'objidx',params.objidx);
    C_total(i)=nansum(abs(Mz).*N0,'all')/nansum(N0(~isnan(Mz)))
end


end

