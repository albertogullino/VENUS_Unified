 Tdu=IdeOon;
  if ick==1
' % Tdu=1; ', keyboard
  end

if length(Li)==0
 Tb=1;
 jsau=0;

else 

 if iTsav==0
  sj=size(Tstor);
  if length(sj)==3
   jsau=sj(3);
  else
   jsau=1;
  end
 else
  jsau=icoustor;
 end

end

 if ifp==-12
  disp(' jsau '), keyboard
 end

if ifp~=-4 & ick==1
' icousav', keyboard
end

%' sono in icrit: chain_vecchio', keyboard
 ficri=find(icrit==1);
 if isfield(Ps,'igacrit')
  if Ps.igacrit==0
  'RESETTO ICRIT'
   ficri=[];
   icrit=icrit*0;
  end
 end
 if dovesono==1
  ficrit=ficri;
 else
  ficrib=ficri;
 end 
 

%ifiez, keyboard

if length(find(icrit>0))>0 & ifiez~=1

%if exist('icrit')==1
 
 %sha_all=shavet(2:end-1); 
 %clear nimi
 %for ks=1:length(nitot)
 % nim=find(nitot(ks,:)>0);
 % nimi(ks)=min(nitot(ks,nim));
 %end 
 %ficr=find(sha_all>3 & nim<2 & aitot>0);
 %jsa=ficr(end);
 %icrit(jsa+1:end)=0;
% ' qui cont nuovo prima', keyboard
if dovesono==1
  jsa=ficri(end); 
  jsau=icoustor; 
else
  jsa=0; 
  jsau=ficri(1)-1; 
end
 % [Gacritn,Trcrit]=Gam_critU(Tstor,Ga1,Mcrit,ficri);
 if ifp==-10
%  ' qui cont nuovo', keyboard
 end
 if dovesono==1
  Ga_old=Ga1;
  Ga_inf=Ga1;
 else
  Ga_old=Ga2;
  Ga_inf=Ga2; 
 end  

%'Gacrit', keyboard

% [Gacrit,Trcrit,Trc,Grc]=Gam_critUn(Tstor,Ga1,Mcrit,ficri,fmlsi);
ficris=ficri;
if igraef_new==2
 fisha=find(shai==6);
 ficri(fisha)=-ficri(fisha);
end 


[Gacrit,Trcrit,Trc,Grc]=Gam_critScattHCG(Ga_inf,Mcrit,ficri,fmlsi,dovesono,icpo);


% [Gacrit,Trcrit]=Gam_critU(Tstor,Ga1,Mcrit,ficri);
 if ifp==-10
%  ' sono in icrit: chain_vecchio DOPO', keyboard
 end
 if dovesono==1
  Trcritu{icpo}=Trcrit;
  Ga1=Gacrit;
 else
  Ga2=Gacrit;
  Trcritb=Trcrit;
 end  
% ' sono in icrit: chain_vecchio', keyboard
else 
 jsa=0;
end

' jsau Mcri ', keyboard
if jsau>0
% jsa=0;
 while jsa<jsau

  jsa=jsa+1;
  nrig=find(icousav==jsa);
  if length(nrig)>1
%   ' nrig in Chain_i ', keyboard
  end
  nstrat=fmlsi(nrig,1);
  if ifp==-11
  disp(' in while 1')
  jsa
  nstrat
  nrig
  pausak
  end
  if ilaymem(nrig)==0
   nstrat=1;
  end
  if nstrat==1
   if iTsav==0
    Tdu=Tstor(:,:,jsa)*Tdu;
   else
%    eval([' load ',nTstof,num2str(jsa)]);
    if ispeed==1
     eval([' load ',Dsav,'\',nTstof,num2str(jsa)]);
    end
    Tdu=Tstof*Tdu;
   end
  else
   nmirro=max(fmlsi(nrig,2))
%   keyboard
   if ifp==-11
    'nmirro'
    nmirro
    pausak
    keyboard
   end
   Tmirro=IdeOon;
   for kmir=1:nstrat
    if iTsav==0
     Tdum=Tstor(:,:,jsa);
    else
     if ispeed==1
      eval([' load ',Dsav,'\',nTstof,num2str(jsa)]);
     end
     Tdum=Tstof;
    end
    Tmirro=Tdum*Tmirro;
    jsa=jsa+1;
    if ifp==-11
     disp(' in while 2')
     jsa
    end
   end

%' nmirro', keyboard
   Pow=Tmirro^nmirro;
   Tdu=Pow*Tdu;
    if ifp==-11
    'Tmirro ', keyboard
     disp(' in while 3')
     jsa
     pausak
    end
   jsa=jsa-1;
  end
%  istr
%  pausak
%  if length(find(istfie==istr))==1
%   icmem=icmem+1;
%   Tmeduf(:,:,icmem)=T;
%   disp(' memorizzo  Tmef')
%   keyboard
%  end
%  jsa
%  pausak
%  if length(find(istfie==jsa))==1
  
%  puf=find(ilaymem==1);
%  if length(puf)==0
%   puf=-1;
%  end
%  if length(istfie)>0
%  if puf==istfie
%   icmem=icmem+1;
%   Tmeduf(:,:,icmem)=Tdu;
%   if ifp~=-4 & ick==1 , disp(' memorizzo  Tmef'), keyboard, end
%  end
 end

else
 if length(Li)>0
 if iTsav==0
  Tdu=Tstor;
 else
  Tdu=Tstof;
 end
 end
end

if ifp==-12
  disp('fine  Tstor '), keyboard
end
%  disp('fine  chain '), keyboard
% clear T Mo
