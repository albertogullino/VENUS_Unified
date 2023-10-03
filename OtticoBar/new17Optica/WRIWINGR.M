NP00=10;
stmin=2.5e-3;

%structure material=gaas alloy=al length=0.05 conc=0.90
%structure material=gaas alloy=al length=0.02 conc=0.90-(0.78/0.02)*d
%structure material=gaas alloy=al length=0.039 conc=0.12
%structure material=gaas alloy=al length=0.02 conc=0.12+(0.78/0.02)*d



fip=find(nomeF=='.');
nomeW=[nomeF(1:fip) 'dev'];

fid = fopen(nomeW,'w');
Ndr=dov;

%' wriwingr', keyboard

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

fid = fopen(nomeW,'w');
fista=find(abs(dov)>0);
icw=fista(1)-1;
icwu=fista(end);
flare=0;

eval(['fprintf(fid,''' cline '\r\n'');']);

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

 if fmol>1 & fmolm~=fmol
% if fmol>1 & fmolm==1
  eval(['fprintf(fid,''' cline '\r\n'');']);
  eval(['fprintf(fid,''' cres '\r\n'');']);
  eval(['fprintf(fid,''' cline '\r\n'');']);
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
  else
   if Ndr(icw)>0
    cN=cNa;
   else
    cN=cNd;
   end
   DOP=abs(Ndr(icw));
  end
%structure material=gaas alloy=al length=0.0066667 conc=0.25

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

fclose(fid);

disp(' ')
disp(' ')
disp( ' wriwingr')
disp(' ')
disp(' ')
%disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ')
%disp(' ')
%disp( ['     Edit and save           ' nomeW])
%disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ')
%disp(' ')
disp(' ')
disp( ' Run SimWindows and save the static EF at the desired voltages '),
disp([' in files with name ' nomeW(1:end-4) '_##, where # are 2 digit = 10 * App. Voltage'])
keyboard
