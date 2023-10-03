function expm_miu(xin,ifx,iacc,Pust);

%tic

global IdeOo IdeOon pMu0u pMu0u1 lKA nk1max nures pMc sim0 ldap
global Oo Odu
%ifx
%' EXp miu in',

Oo=xin(Pust,Pust);

 if ifx==-1
  Oo=IdeOon+Oo;
 elseif ifx==-10
%   eM=expm_3(x);
   expm_3u;
 elseif ifx==0
% ' ifx ', keyboard
  if length(iacc)<=1
   expm_3u;
  else
   Oo=IdeOon;

%   Kplot1=Oo;
%   an_mate
   ' an_mate ', keyboard

   sia=size(iacc);
   for ki=1:sia(1)
    puu=[];
    fi=find(iacc(ki,:)~=0);
    iaccv=iacc(ki,fi);
    if length(iaccv)>0
     for ipu=iaccv
      for ipui=iaccv
       add=((ipu-1)*2*lKA+ipui-1)*nk1max;
       puu=[puu pMc+add];
      end
     end
     puu=sort(puu);
     sim=sim0*length(iaccv);
     xx=reshape(Oo(puu),sim,sim);
     Oo(puu)=expm3(xx);
%     mapab(log10(x)), pausak
%     mapab(log10(xx)), pausak
%     mapab((eM)), pausak
    end
   end
  end

 elseif ifx==1

%  sx=size(x);
%  l=sx(1)/2;
%  x2=diag(x,l);
%  xd=diag(x,0);
%  x1=xd(1:l);
%  x4=xd(l+1:2*l);
%  x3=diag(x,-l);
%  for k=1:l
%   xx=[x1(k) x2(k); x3(k) x4(k)];
%   em(k,:)=reshape(expm(xx),1,4);
%  end
%  di=[em(:,1); em(:,4)];
%  x2=em(:,3);
%  x3=em(:,2);
%  eM=diag(di)+diag(x2,l)+diag(x3,-l);

  Odu=IdeOo;
  for ki=1:lKA
   puu=pMu0u1+ki+(ki-1)*lKA*2;
   xx=reshape(xin(puu),2,2);
   xm=expm3(xx);
   Odu(puu)=xm;
  end
  Oo=Odu(Pust,Pust);

 elseif ifx==2
  Odu=IdeOo;
  for ki=1:nk1max
   puu=pMu0u+ki+(ki-1)*2*lKA;
   xx=reshape(xin(puu),nures,nures);
   xm=expm3(xx);
   Odu(puu)=xm;
  end


  Oo=Odu(Pust,Pust);

 end


%' EXp miu out',
%toc,
%pausak
