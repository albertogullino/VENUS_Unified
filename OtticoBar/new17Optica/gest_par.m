function [cPAR,cPver]=gest_par(nomestr,ifp)

%'gest_par i', keyboard
if nargin==1
    ifp=-10;
end
% ip=findstr('.str',nomestr_);
% if length(ip)==0
%  disp' errore gest_par')
% else
  nopar=[nomestr,'.par'];
% end
% nopamat=[nomestr,'par'];
if ifp~=-4
  %  eval(['type ',nopar]); % par visualized in Command window
%     disp('')
%     disp('')
%     disp('   &&&&&&&&&&&&&&&&&&& ')
%     disp('')
%     disp('')
%     ichange=input(' 1 per cambiare i parametri, ENTER per proseguire ');
%     disp('')
%     disp('')
ichange=[];
else
    ichange=[];
end
if length(ichange)==1
    disp('')
    disp('')
    dis=[' EDIT, MODIFICA E SALVA          ------->        ',nopar];
    disp(dis)
    dis=[' NB: solo i parametri con il carattere @ sono considerati '];
    disp(dis)
    if ifp>-4
        keyboard
    end
end
fid=fopen(nopar,'r');
ico=1;
while feof(fid)==0
    Nline=fgetl(fid);
    fc=findstr(Nline,'@');
    if length(fc)==1;
        [p1,remainder]=strtok(Nline);
        p2=Nline(fc+2:end);
        %       disp(Nline)
        cPver(ico,1)=str2num(p1);
        cPdu(ico,1:length(p2)+1)=[p2,' '];
        Nadd=['          ',num2str(ico),'  ',Nline];
        cDisp(ico,1:length(Nadd))=Nadd;
        ico=ico+1;
    end
end
if exist('cPdu')
    s=size(cPdu);
    scol=s(2);
    cPAR=reshape(cPdu',prod(s),1)';
else
    cPAR=[];
    cPver=[];
end

if ifp~=-4
%     [' PARV (active parameter marked by @) must have the same order as the listed parameters ']
%     [' Index PARV,      Par,      Layer,  Tran_sect,   Par. type,       shape,       short_name']
%     disp( cDisp)
    
%     pausak
end

%'gest_par u', keyboard
