NP00=10;
stmin=2.5e-3;



fip=find(nomeF=='.');
%nomeW=[nomeF(1:fip) 'dev'];
nomeUP=[nomeF(1:fip-1) 'UP.dev'];
nomeDW=[nomeF(1:fip-1) 'DW.dev'];
nomedu=[nomeF(1:fip-1) 'du.dev'];

fidUP = fopen(nomeUP,'w');
fidDW = fopen(nomeDW,'w');
fiddu = fopen(nomedu,'w');
Ndr=dov;

%' wriwimir', keyboard

LW=dvw*1e-3;
cx=xmw(:,1);
fic=find(cx<1e-6);
cx(fic)=0;


cst='structure material=gaas alloy=al length=';
cle=' conc='                                  ;
cgr='grid length=';
cpt=' points=';
cdp='doping length=';
cNa=' Na=';
cNd=' Nd=';
cres='repeat start';
cre='repeat=';
crebu='region bulk length=';
creqw='region qw length=';
cline='                                        ';

fista=find(abs(dov)>0);
icw=fista(1)-1;
icwu=fista(end);
flare=0;

eval(['fprintf(fidUP,''' cline '\r\n'');']);
eval(['fprintf(fidDW,''' cline '\r\n'');']);
eval(['fprintf(fiddu,''' cline '\r\n'');']);

flag_up=0;
flag_dw=0;
flag_du=0;

fid=fidUP;

%while icw<icwu-1
while icw<icwu
 icw=icw+1;
if LW(icw)>0

 fmol=abs(fsw(icw));

 if icw>1
  fmolm=abs(fsw(icw-1));
 else
  fmolm=1;
 end

 if icw<length(LW)
  fmolp=abs(fsw(icw+1));
 else
  fmolp=1;
 end



 thick=LW(icw);
 NP0=max([fix(thick/stmin) NP00]);


  if Ndr(icw)==0
   if Ndr(icw+1)>0
    cN=cNa;
   else
    cN=cNd;
   end
   DOP=abs(Ndr(icw+1));
   DOPs=(Ndr(icw+1));
  else
   if Ndr(icw)>0
    cN=cNa;
   else
    cN=cNd;
   end
   DOP=abs(Ndr(icw));
   DOPs=(Ndr(icw));
  end

  if DOPs~=0
   if DOPs>0
    flag_up=1;
   end
   if DOPs<0
    flag_dw=1;
    flag_up=0;
   end
  end
  if flag_up==1 & abs(DOPs)<=1e17
   fid=fiddu;
  end
  if flag_dw==1
   fid=fidDW;
  end

 if fmol>1 & fmolm~=fmol
  eval(['fprintf(fid,''' cline '\r\n'');']);
  eval(['fprintf(fid,''' cres '\r\n'');']);
  eval(['fprintf(fid,''' cline '\r\n'');']);
 end

   eval(['fprintf(fid,''' cline '\r\n'');']);
  if imag(cx(icw))==0
   eval(['fprintf(fid,''' cst '%0.5g' cle '%0.5g\r\n'',LW(icw),cx(icw));']);
  else
   cx1=real(cx(icw));
   cx2=imag(cx(icw));
   Ccx1=num2str(cx1);
   Ccoef=num2str(abs((cx2-cx1)/LW(icw)));
   if cx1>cx2
    sig='-';
   else
    sig='+';
   end
   stru=[cst,num2str(LW(icw)),cle,Ccx1,sig,Ccoef,'*d \r\n'];
%   'stru'
%   keyboard

   eval(['fprintf(fid,stru);']);
  end
   eval(['fprintf(fid,''' cgr '%0.5g' cpt '%0.5g\r\n'',LW(icw),NP0);']);
   eval(['fprintf(fid,''' cdp '%0.5g' cN '%0.5g\r\n'',LW(icw),DOP);']);
   if fsw(icw)==-1
    eval(['fprintf(fid,''' creqw '%0.5g\r\n'',LW(icw));']);
   else
    eval(['fprintf(fid,''' crebu '%0.5g\r\n'',LW(icw));']);
   end
   eval(['fprintf(fid,''' cline '\r\n'');']);
% [fmol fmolp]

% if fmolp==1 & fmol>1
 if fmolp~=fmol & fmol>1
  npair=fmol;
  eval(['fprintf(fid,''' cline '\r\n'');']);
  eval(['fprintf(fid,''' cre '%0.5g\r\n'',npair);']);
  eval(['fprintf(fid,''' cline '\r\n'');']);
% 'fmol'
% [fmol fmolp]
% pausak
 end

end
end

fclose(fiddu);
fclose(fidDW);
fclose(fidUP);

disp(' ')
disp(' ')
disp( ' wriwimir: this is only for the resistance of the mirrors')
disp(' ')
disp(' ')
%disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ')
%disp(' ')
%disp( ['     Edit and save           ' nomeW])
%disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ')
%disp(' ')
disp(' ')
%disp( ' Run SimWindows and save the static EF at the desired voltages '),
%disp([' in files with name ' nomeW(1:end-4) '_##, where # are 2 digit = 10 * App. Voltage'])
%keyboard
