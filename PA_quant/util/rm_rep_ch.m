function str=rm_rep_ch(str,ch)

old_str=str;
while str
    
    str=strrep(str,[ch,ch],ch);
    if strcmp(str,old_str)
        break;
    else
        old_str=str;
    end
    

end


