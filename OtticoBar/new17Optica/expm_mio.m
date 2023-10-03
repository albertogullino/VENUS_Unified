function [x,ier]=expm_mio(xin,ifx,iacc,Pust,PUrid);

x=xin(Pust,Pust);

global IdeOo IdeOon pMu0u pMu0u1 lKA nk1max nures pMc sim0 ldap ive
global ldapu MPUa
global ired_ret

%if length(iacc)>0
%
%' EXp mio iama: ifx',ifx, keyboard
%end

%' EXp mio in: ifx',ifx, pausak
  [V,autov]=eig(x);
  au=diag(autov);
  maau=max(real(au));
  %pausak
  if maau>.1
   ifx=-1000;
  end

 if ifx==-1000
   expm_3p;  
 elseif ifx==-1 
  x=IdeOon+x;
 elseif ifx==-10
%   eM=expm_3(x);
   expm_3;
 elseif ifx==0
  if length(iacc)<=1
   expm_3;
  else
   xex=x;
   x=IdeOon;

%   an_mate
%   ' sono qui per anmate ',
   idisp=0;
   if idisp==1

   PU=zeros(168,168);
   for ic=1:sia
    puu=PUrid{ic};
    sim=sqrt(length(puu));
%    xx=reshape(xex(puu),sim,sim);
    xx=PU;
    xx(puu)=1e5;
    mapab(xx), pausak
%    x(puu)=expm3(xx);
   end
   end



   sia=length(PUrid);
   for ic=1:sia
    puu=PUrid{ic};
    sim=sqrt(length(puu));
%    xxdu=xex;
%    xxdu(puu)=1e5;
%    mapab(xxdu), pausak
%    xxdu=xex;
%    puu1=PUrid{2};
%    xxdu(puu1)=0;
%    mapab(xxdu),   pausak
%    ' qui planare ', keyboard
     io=1;
%     'tscatt2', keyboard
     xx=reshape(xex(puu),sim,sim);
     if ired_ret~=1
      x(puu)=expm(xx);
     else 
%      [xve,Dau] = eig(xx);
%      xax = xve * diag(exp(diag(Dau))) / xve;
%'tscatt2', keyboard
      Tscatt0
%      Tscatt1
%      Tscatt
      x(puu)=xax;
%' in exp_mio ', keyboard
     end    
   end
   iversem=0;
   if iversem==1
    xver=expm3(xin);
%    xver=expm(xin);
    mapab(xver-x), pausak
   end

%   for ki=1:sia(1)
%    puu=[];
%    fi=find(iacc(ki,:)~=0);
%    iaccv=iacc(ki,fi);
%    if length(iaccv)>0
%     for ipu=iaccv
%      for ipui=iaccv
%       add=((ipu-1)*2*lKA+ipui-1)*nk1max;
%       puu=[puu pMc+add];
%      end
%     end
%     puu=sort(puu);
%     sim=sim0*length(iaccv);
%     xx=reshape(x(puu),sim,sim);
%     x(puu)=expm3(xx);
%%     mapab(log10(x)), pausak
%%     mapab(log10(xx)), pausak
%%     mapab((eM)), pausak
%    end
%   end


  end

 elseif ifx==1
%%%%%%%qui  'ifx=1 ', keyboard

  x=IdeOo;
  for ki=1:lKA
   puu=pMu0u1+ki+(ki-1)*lKA*2;
   xx=reshape(xin(puu),2,2);
%   xm=expm3(xx);
%   xm=expm(xx);
        if ired_ret~=1
%         xx=reshape(xex(puu),2,2);
         xm=expm(xx);
        else 
   %      [xve,Dau] = eig(xx);
   %      xax = xve * diag(exp(diag(Dau))) / xve;
%   'tscatt2', keyboard
%         Tscatt2
%         Tscatt1
%         Tscatt
         Tscatt0
         xm=xax;
   %' in exp_mio ', keyboard
     end    
   
   x(puu)=xm;
%   ki
%   pausak

  end
  x=x(Pust,Pust);

 elseif ifx==2

%  ' EXp mio in: ifx',ifx, keyboard
  x=IdeOo;
  ipri=1;
  isal=1;
  ipri=0;
  isal=0;
  if isal==1
   ifx
  end
 if ive==1

  for ki=1:nk1max
   puu=pMu0u+ki+(ki-1)*2*lKA;
   xx=reshape(xin(puu),nures,nures);
%   xm=expm3(xx);
   xm=expm(xx);
   x(puu)=xm;
   if isal==1
    xp=xm;
    save xv xp
    ki, pausak
   end
   if ipri==1
    px=xin;
    px(puu)=px(puu)*100;
    ki
    map(px), pausak
   end
  end
   if isal==1
    xp=x;
    save xv xp
    'dopo', pausak
   end
 else
  for ki=1:nk1max
   puu=pMu0u{ki};
   s=length(puu);
   nur=sqrt(s);
   xx=reshape(xin(puu),nur,nur);
%   xm=expm3(xx);
   xm=expm(xx);
   x(puu)=xm;
   if isal==1
    xp=xm;
    save xn xp
    ki, pausak
   end
   if ipri==1
    px=xin;
    px(puu)=px(puu)*100;
    ki
    map(px), pausak
   end
  end
   if isal==1
    xp=x;
    save xn xp
    'dopo', pausak
   end
 end


  x=x(Pust,Pust);

%  ' dopo x ', pausak

 end

DE=det(x);
%pausak
ier=0;
if abs(DE-1)>1e-10
 ier=1;
end

%' EXp mio out'
%toc,
%pausak

%ifx
%' EXp mio out: ifx',ifx, pausak
