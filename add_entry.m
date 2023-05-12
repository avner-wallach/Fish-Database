function [db] = add_entry(db,base)
varnames=db.Properties.VariableNames(end);
dlgtitle = ['New Entry'];
prompt = varnames;
dims = [1 35];
definput = {''};
answer = inputdlg(prompt,dlgtitle,dims,definput);
db=[db;base str2num(answer{1})];
end

