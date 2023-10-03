iplover=0;
PUridp=[];
PUridm=[];
%global PUridm PUridp
iaccp=[];
iaccm=[];
%return
if ~exist('farm')==1
 farm=1;
end
if farm>1
 return
end
%if strcmp(Ev_or_Od,'Odd')==1
% return
%end
%if sha~=1 | (sha==1 & is_even(mm)==1) | iany==1 | ianys==1
%if sha~=1 | iany==1 | ianys==1
%if sha~=1
% return
%end
%' passo ??? ', keyboard

if iLP==1
 npias=numodi;
else
 npias=2*numodi;
end
lKAp=numodi*npias;
clear Pias1 Pias2
ld1=[1:nk1max];
M1=max(max(abs(Kplot1)));
M2=max(max(abs(Kplot2)));
if M1==0
 M1=1e-6;
end
if M2==0
 M2=1e-6;
end

%  if iLP==0
%   nktot=2*nk1max*numodi;
%  else
%   nktot=1*nk1max*numodi;
%  end
  nktot=Pust(end);
%  nktot=length(Kplot1);
  VPU=1:nktot^2;
  MPU=reshape(VPU,nktot,nktot);
%  'MPU', keyboard

for ir=1:npias
%  pMp=[];
%  for kdu=ld1
%   pMp=[pMp ld1+(ir-1)*nk1max+(kdu-1)*lKAp];
%  end
 for ic=1:npias
%  pMpd=pMp+(ic-1)*nk1max*lKAp;
  Mdu=MPU;
  Mdu(ldapu(ir)+1:ldapu(ir+1),ldapu(ic)+1:ldapu(ic+1))=NaN;
  pMpd=find(isnan(Mdu)==1);
%
%  lP=length(PUridm);
%  for kpl=1:lP
%   Kdu=Kplot1;
%   pMpd=PUridm{kpl};
%   Kdu(pMpd)=1e5;
%   map(Kdu), %pausak
%  end
%
%
%
%
  m=mean(abs(Kplot1(pMpd)));
  if iplover>1
   Kdu=Kplot1;
   Kdu(pMpd)=1e5;
   map(Kdu),
   m/M1,  pausak
  end
  if m/M1<1e-6
   Pias1(ir,ic)=0;
  else
   Pias1(ir,ic)=1;
  end
  m=mean(abs(Kplot2(pMpd)));
%  m/M2
%  pausak
  if m/M2<1e-6
   Pias2(ir,ic)=0;
  else
   Pias2(ir,ic)=1;
  end
  if ic==ir
   Pias1(ir,ic)=1;
   Pias2(ir,ic)=1;
  end
 end
end

%' fine PIAS', keyboard

clear iaccp iaccm
clear iaccpd iaccmd


for ir=1:npias
 iacc=find(Pias1(ir,:)==1);
 if length(iacc)~=npias & length(iacc)~=0
  iaccpd(ir,1:length(iacc))=iacc;
 end
end

for ir=1:npias
 iacc=find(Pias2(ir,:)==1);
 if length(iacc)~=npias & length(iacc)~=0
  iaccmd(ir,1:length(iacc))=iacc;
 end
end

if exist('iaccpd')
 [su,iso]=sort(iaccpd(:,1));
 iaccpd=iaccpd(iso,:);
 fiv=find([1; diff(iaccpd(:,1))]~=0);
 iaccpdu=iaccpd(fiv,:);
 sa=size(iaccpdu);
 iaccp=[];
 for ka=1:sa(1)
  if length(find(iaccpdu(ka,:)~=0))>0
   iaccp=[iaccp; iaccpdu(ka,:)];
  end
 end
else
 iaccp=[];
end

if exist('iaccmd')
 iaccpd=iaccmd;
 [su,iso]=sort(iaccpd(:,1));
 iaccpd=iaccpd(iso,:);
 fiv=find([1; diff(iaccpd(:,1))]~=0);
 iaccpdu=iaccpd(fiv,:);
 sa=size(iaccpdu);
 iaccm=[];
 for ka=1:sa(1)
  if length(find(iaccpdu(ka,:)~=0))>0
   iaccm=[iaccm; iaccpdu(ka,:)];
  end
 end
else
 iaccm=[];
end

%if length(iaccp)>0
% fi0=find(iaccp==0);
% if length(fi0)>0
%  pu=fi0-length(iaccp);
%  iaccp(fi0)=iaccp(pu);
% end
%end
%
%if length(iaccm)>0
% fi0=find(iaccm==0);
% if length(fi0)>0
%  pu=fi0-length(iaccm);
%  iaccm(fi0)=iaccm(pu);
% end
%end

%' qui ', keyboard

%  if iLP==0
%   nktot=4*nk1max*numodi;
%  else
%   nktot=2*nk1max*numodi;
%  end

% if iLP==0
  nktot=2*length(Kplot1);
% else
%  nktot=length(Kplot1);
% end
  nktot=Pust(end);
  VPU=1:nktot^2;
  MPU=reshape(VPU,nktot,nktot);

pros=prod(size(MPU));
ldapuu=[ldapu ldapu(2:end)+ldapu(end)];

if length(iaccm)==0
 PUridm=[];
else

 for irid=1:length(iaccm(:,1))
 iadu=iaccm(irid,:);
 iadu1=iadu(find(iadu>0));
 iaccdu=[iadu1 iadu1+npias];
% iaccdu=[iaccm(irid,:) iaccm(irid,:)+npias];
% pausak
 pMdu=[];
 for ir=iaccdu
  for ic=iaccdu
   Mdu=MPU;
   if ic>0 & ir>0
    Mdu(ldapuu(ir)+1:ldapuu(ir+1),ldapuu(ic)+1:ldapuu(ic+1))=NaN;
%    map(Mdu), pausak
    pMpd=find(isnan(Mdu)==1);
    pMdu=[pMdu; pMpd];
   end
  end
 end
  if iplover==1
   Mdu=MPU;
   Mdu(pMdu)=NaN;
   map(Mdu), pausak
  end
%  PUridm{irid}=pMdu;
  PUridm{irid}=sort(pMdu);
 end
end

if length(iaccp)==0
 PUridp=[];
else

 for irid=1:length(iaccp(:,1))
% iaccdu=[iaccp(irid,:) iaccp(irid,:)+npias];
 iadu=iaccp(irid,:);
% pausak
 iadu1=iadu(find(iadu>0));
 iaccdu=[iadu1 iadu1+npias];
 pMdu=[];
 for ir=iaccdu
  for ic=iaccdu
   if ic>0 & ir>0
    Mdu=MPU;
    Mdu(ldapuu(ir)+1:ldapuu(ir+1),ldapuu(ic)+1:ldapuu(ic+1))=NaN;
    pMpd=find(isnan(Mdu)==1);
    pMdu=[pMdu; pMpd];
   end
  end
 end
%  PUridp{irid}=pMdu;
  PUridp{irid}=sort(pMdu);
 end
end

%'an_mat '
% keyboard

if iplover==1
 xex=[Kplot1 Kplot1; Kplot1 Kplot1];
 mapab(xex);
   PUrid=PUridp;
   sia=length(PUrid);
   for ic=1:sia
    puu=PUrid{ic};
    sim=sqrt(length(puu));
    xdis=xex;
    xdis(puu)=xdis(puu)*50;
    mapab(xdis); pausak
   end

 ' an_mat ', keyboard
end
