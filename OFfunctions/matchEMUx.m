function hitEMU = matchEMUx(EMUname, EMUtag, EMUlist)
if isempty(EMUlist)
    hitEMU = [];
    return
end

hitEMU = strcmp(EMUname,EMUlist(:,1));
for i = find(hitEMU)'
    if ~all(EMUlist{i,2}==EMUtag)
        hitEMU(i) = false;
    end
end