function [ net,rsq,sig ] = ni_network_new(inputs,targets,mode,dims,nhidden)
%NI_NETWROK generate and train negative image regression. hidden layer is
%random and static, output layer is performing gradient decent
if(nargin<5)
    if(strcmp(mode,'BP'))
        nhidden=100;
    else
        nhidden=1e3;
    end
end
est=0;
outlayer='tansig';
% outlayer='purelin';

if(nargin<4 | ~numel(dims))
    dims=size(inputs,1);
end
%% generate backprop network
if(strcmp(mode,'BP'))
    net = feedforwardnet(nhidden);
    net = fitnet(nhidden);
%     net.layers{2}.transferFcn = outlayer;
    % train network
    net = train(net,inputs,targets,'useParallel','no','useGPU','no');
    outputs=net(inputs);
    rsq=nancorr(outputs(:),targets(:))^2;
    sig=nanstd(outputs(:)-targets(:));
    return;
end
%% generate granular stage
maxD=inf; %maximum number of inputs to GC
y=[]; x=[];
k=0; I=0;
for d=1:numel(dims)    
    p=dims(d)/sum(dims);
    for i=1:round(nhidden*p) %each row = one GC
        J=randperm(dims(d),randi(min(dims(d),maxD)))+k;
        y=[y;J'];
        x=[x;(I+i*ones(size(J')))];
    end
    k=k+dims(d);
    I=I+i;
end
w=randn(size(x));
W=sparse(x,y,w);    
B=randn(nhidden,1);
mf=W*inputs+B; %'synaptic input' into GC
if(nhidden<=2e3)
    gc=tansig(mf);  %GC output activity
else
    gc=zeros(size(mf));
    for i=1:size(mf,2)
        gc(:,i)=tansig(mf(:,i));
    end
end
%% generate network[
net1 = fitnet(1,'trainbr');
% net.layers{1}.transferFcn = outlayer;
net1.inputs{1}.processFcns={};
net1.outputs{2}.processFcns={'mapminmax'};
net1=configure(net1,gc,targets);
net1.biases{2}.learn=0;
net1.b{2}=0;
net1.layerWeights{2,1}.learn=0;
net1.LW{2,1}=1;
%% train network
net1 = train(net1,gc,targets,'useParallel','no');%,'useGPU','yes');

%% generate output network
net = fitnet(nhidden);
net.inputs{1}.processFcns={};
net.outputs{2}.processFcns={'mapminmax'};
net.layers{2}.transferFcn = outlayer;
net = configure(net,inputs,targets);
net.iw{1,1}=W;
net.b{1}=B;
net.b{2}=net1.b{1};
net.LW{2,1}=net1.iw{1,1};

outputs=net(inputs);
rsq=nancorr(outputs(:),targets(:))^2;
sig=nanstd(outputs(:)-targets(:));

%% estimate degrees of freedom
if(est)
%     enet = feedforwardnet(nhidden);
%     enet.layers{2}.transferFcn = outlayer;
    N=size(inputs,1);    
    for i=1:N %go over each input
        I=inputs;
        I(i,:)=nanmedian(I(i,:));
%         I(i,:)=I(i,randperm(size(I,2))); %shuffle
        outputs=net(I);
        prsq(i)=(rsq-nancorr(outputs(:),targets(:))^2)/rsq;
%         ind=setdiff(1:N,i);
%         enet = train(enet,inputs(ind,:),targets,'useParallel','yes');%,'useGPU','yes');
%         eoutputs=enet(inputs(ind,:));
%         ersq=nancorr(eoutputs(:),targets(:))^2;
%         prsq(i)=(rsq-ersq)/rsq;
    end    
else
    prsq=[];    
end

end