function [ group ] = process_group_models(eod0,file0,ops,group)
%PROCESS_GROUP_MODELS get all stats for group

setenv('PXLSIZE','0.423'); %pxl to mm conversion

groupnum=group.groupname;

process;


    function process
        if(~numel(segnum))
            group.model=[];
            return;
        end

        eod=eod0(segnum);
        file=file0(segnum);

        %get cols
        lfpcol=find(cellfun(@(x) numel(strfind(x,['lfp',num2str(groupnum)])),eod.fnames));
        if(~isfield(group,'model'))
            group.model=get_regression_models(eod,file,lfpcol,'objidx',objidx);
        end
        [ group.decoding ] = evaluate_regression_models(eod,file,group.model,'plotting',0,'objidx',objidx);
    end    

end
