function P=prodmat(xin,yin,ifx,Pust);
global IdeOo pMu0u pMu0u1 lKA nk1max nures pMc ive
%' prodmat'
%tic
global PUrid
%ifx=0;

if length(xin)>1000
 P=xin*yin;
 return
end

x=IdeOo;
y=IdeOo;
x(Pust,Pust)=xin;
y(Pust,Pust)=yin;


six=size(xin);
siy=size(yin);
iacc=PUrid;
%iacc=0;
tic
if six~=siy
 if prod(six)==1 | prod(siy)==1
  P=xin*yin;
 else
  disp(' prodmat.m:  matrices of different dim ')
  keyboard
 end
 return
else
 if prod(six)==1
  P=xin*yin;
  return
 end
end

 if ifx==0

  if length(iacc)<=1
   P=xin*yin;
  else
   P=zeros(size(xin));
    for ki=1:length(PUrid)
     puu=PUrid{ki};
     sim=sqrt(length(puu));
     xx=reshape(x(puu),sim,sim);
     yy=reshape(y(puu),sim,sim);
     xm=xx*yy;
     P(puu)=xm;
    end
   end

 elseif ifx==1


  P=IdeOo;
  for ki=1:lKA
   puu=pMu0u1+ki+(ki-1)*lKA*2;
   xx=reshape(x(puu),2,2);
   yy=reshape(y(puu),2,2);
   xm=xx*yy;
   P(puu)=xm;
  end
   P=P(Pust,Pust);

 elseif ifx>=2

  isal=0;
  P=IdeOo;
 if ive==1
   if isal==1
   xp=x;
   yp=y;
   save xv xp yp
   ' prima prod vec ', pausak
   end
  for ki=1:nk1max
   puu=pMu0u+ki+(ki-1)*2*lKA;
   xx=reshape(x(puu),nures,nures);
   yy=reshape(y(puu),nures,nures);
   xm=xx*yy;
   P(puu)=xm;
  end
   if isal==1
   xp=P(Pust,Pust);
   save xv xp yp
   ' dopo prod vec ', pausak
   end
 else
   if isal==1
   xp=x;
   yp=y;
   save xn xp
   ' prima prod nuo ', pausak
   end
  for ki=1:nk1max
   puu=pMu0u{ki};
   s=length(puu);
   nur=sqrt(s);
   xx=reshape(x(puu),nur,nur);
   yy=reshape(y(puu),nur,nur);
   xm=xx*yy;
   P(puu)=xm;
  end
   if isal==1
   xp=P(Pust,Pust);
   save xn xp
   ' dopo prod nuo ', pausak
   end
 end


   P=P(Pust,Pust);

 end
if ifx==100
 toc
'prodmat'
keyboard
end

%toc,
%pausak
