NP00=10;
stmin=2.5e-3;



fip=find(nomeF=='.');
nomeW=[nomeF(1:fip) 'dev'];


fid=find(diff(dov)~=0)+1;
fid=fid(2:end);
fi4=fid;
Ndr=dov;



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
icw=1;
flarepe=0;
subu=0;
subug=0;
subud=0;
subuto=0;

while icw<length(LW)-1
 icw=icw+1;
 fmol1=1;
 fmol=abs(fsw(icw));
 fmol2=abs(fsw(icw));
 fst0=fsw(icw-1);
 fst1=fsw(icw);
 'wre', pausak
 if ( (fst0>1 & fst1>1) & fst0~=fst1)
  fla2re=1;
 else
  fla2re=0;
 end
 if length(find(fi4-(icw)==0))>0

   if Ndr(icw-1)>0
    cN=cNa;
   else
    cN=cNd;
   end

   eval(['fprintf(fid,''' cline '\n'');']);
%   disp(' flarepe, subud'),
%   [flarepe, subud]
%   keyboard

   if flarepe==1
    fdiv=fmol;
   else
    fdiv=1;
   end
   eval(['fprintf(fid,''' cdp '%0.5g' cN '%0.5g\n'',subud/fdiv,abs(Ndr(icw-1)));']);
   eval(['fprintf(fid,''' cline '\n'');']);
   subud=0;
 end

 if abs(fsw(icw))==1 & flarepe==0
  if fsw(icw)==-1
   eval(['fprintf(fid,''' cline '\n'');']);
   NP0=max([fix(subug/stmin) NP00]);
   eval(['fprintf(fid,''' cgr '%0.5g' cpt '%0.5g\n'',subug,NP0);']);
   fmol=0;
   fmol1=0;

   if Ndr(icw)>0
    cN=cNa;
   else
    cN=cNd;
   end

   subu=0;
   subug=0;

   eval(['fprintf(fid,''' cline '\n'');']);
   eval(['fprintf(fid,''' crebu '%0.5g\n'',subuto);']);
   subuto=0;
   eval(['fprintf(fid,''' creqw '%0.5g\n'',LW(icw));']);

  end

  eval(['fprintf(fid,''' cst '%0.5g' cle '%0.5g\n'',LW(icw),cx(icw));']);


 elseif fsw(icw)==1 & flarepe==1
  flarepe=0;

  eval(['fprintf(fid,''' cline '\n'');']);
  NP0=max([fix(subug/stmin) NP00]);
  eval(['fprintf(fid,''' cgr '%0.5g' cpt '%0.5g\n'',subug,NP0*2);']);

  subu=0;
  subug=0;
  eval(['fprintf(fid,''' cre '%0.5g\n'',npair);']);
  eval(['fprintf(fid,''' cline '\n'');']);
  eval(['fprintf(fid,''' cst '%0.5g' cle '%0.5g\n'',LW(icw),cx(icw));']);

 elseif fsw(icw)~=1 & flarepe==0
  npair=fsw(icw);
  eval(['fprintf(fid,''' cline '\n'');']);
   if Ndr(icw-1)>0
    cN=cNa;
   else
    cN=cNd;
   end
   NP0=max([fix(subug/stmin) NP00]);
   eval(['fprintf(fid,''' cgr '%0.5g' cpt '%0.5g\n'',subug,NP0);']);
  subu=0;
  subug=0;
  eval(['fprintf(fid,''' cline '\n'');']);
  eval(['fprintf(fid,''' cres '\n'');']);
  eval(['fprintf(fid,''' cst '%0.5g' cle '%0.5g\n'',LW(icw),cx(icw));']);
  flarepe=1;

 elseif fsw(icw)~=1 & flarepe==1
  if fla2re==1
   npair=fsw(icw-1);
   NP0=max([fix(subug/stmin) NP00]);
   eval(['fprintf(fid,''' cgr '%0.5g' cpt '%0.5g\n'',subug,NP0);']);
   subu=0;
   subug=0;
   eval(['fprintf(fid,''' cre '%0.5g\n'',npair);']);
   eval(['fprintf(fid,''' cline '\n'');']);
   eval(['fprintf(fid,''' cres '\n'');']);
   npair=fsw(icw);
  end
  eval(['fprintf(fid,''' cst '%0.5g' cle '%0.5g\n'',LW(icw),cx(icw));']);

 end

 su=LW(icw)*fmol1;
 sug=LW(icw);
 suo=LW(icw)*fmol;
 sud=LW(icw)*fmol2;
 subu=subu+su;
 subug=subug+sug;
 subud=subud+sud;
 subuto=subuto+suo;

end

   eval(['fprintf(fid,''' cline '\n'');']);
   NP0=max([fix(subug/stmin) NP00]);
   eval(['fprintf(fid,''' cgr '%0.5g' cpt '%0.5g\n'',subug,NP0);']);

   if Ndr(icw)>0
    cN=cNa;
   else
    cN=cNd;
   end

   eval(['fprintf(fid,''' cdp '%0.5g' cN '%0.5g\n'',subud,abs(Ndr(icw)));']);
   subu=0;
   subug=0;

   eval(['fprintf(fid,''' cline '\n'');']);
%   eval(['fprintf(fid,''' crebu '%0.5g\n'',subuto);']);
   eval(['fprintf(fid,''' crebu '%0.5g'',subuto);']);


fclose(fid);

disp( ' wriwinu')
disp(' ')
disp(' ')
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ')
disp(' ')
disp( ['     Edit and save           ' nomeW])
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ')
disp(' ')
disp(' ')
disp( ' then run SimWindows and save the static EF at the desired voltages '),
disp([' in files with name ' nomeW(1:end-4) '_##, where # are 2 digit = 10 * App. Voltage'])
keyboard
