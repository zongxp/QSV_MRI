function res=filename(a)


[res,name,suf]=fileparts(a);

if isempty(res)
    res=a;
else
    res=[name,suf];
end

    