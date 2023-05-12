function [] = plot_vars(group,varnames)
if(~iscell(varnames))
    varnames={varnames};
end
N=numel(varnames);
C=ceil(sqrt(numel(varnames)));
R=ceil(N/C);

figure;
set(gcf,'WindowState','maximized');
[ha, pos] = tight_subplot(R, C, [.1 .1],[.1 .1],[.1 .1]);
for n=1:N
    axes(ha(n));
    if(strfind(varnames{n},'Lfp'))
        bar([cell2mat({group.units.(varnames{n})})]);
        ha(n).XTickLabel={group.units.uname};
    else
        bar([group.lfp.(varnames{n}),cell2mat({group.units.(varnames{n})})]);
        ha(n).XTickLabel=['LFP' {group.units.uname}];
    end
    title(varnames{n},'Interpreter','none')
end
end

