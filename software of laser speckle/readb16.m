function I=b16reader(filename)
[fid,message]=fopen(filename,'rb');
if(fid==-1)
    disp(message)
    return
end

fseek(fid,512+512,'bof');
i=fread(fid,[696 520],'uint16=>uint16');
I=i';
fclose(fid);
