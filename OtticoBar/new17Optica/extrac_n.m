function n=extract_n(Nline,chs)
%         chs='diam=';
         ise=length(chs)+1;
         ipu=findstr(Nline,chs);
        [chr,remainder]=strtok(Nline(ipu:end));
   if strcmp(chs,'Mateq=')==1
    n=chr(ise:end);
   else 
    n=str2num(chr(ise:end));
   end