function [experiment,unit] = generate_expstructs()
% construct for database
unit=struct('uname',0,'baseline',nan,'spmode',nan);
group=struct('groupname',0,'wallseg',0,'brassseg',0);
location=struct('locname',nan,'groups',group);
experiment=struct('expname','','locations',location);
end

