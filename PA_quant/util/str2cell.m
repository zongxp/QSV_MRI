function str = str2cell(str)

if ~iscell(str)
    temp = str;
    str = cell(1,1);
    str{1} = temp;
end