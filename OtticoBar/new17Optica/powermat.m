function P=powermat(xin,p,ifx,Pust);

global IdeOo pMu0u pMu0u1 lKA nk1max nures pMc ive
global PUrid
%ifx=0;

%  P=xin^p;
%  return
%' power '
%tic
x=IdeOo;
y=IdeOo;
x(Pust,Pust)=xin;


iacc=PUrid;
%iacc=0;
six=size(xin);
 if prod(six)==1
  P=xin^p;
  return
 end

 if ifx==0
  if length(iacc)<=1
   P=xin^p;
  else
%   P=IdeOo;
%   sia=size(iacc);
%   for ki=sia(1)
%    puu=[];
%    for ipu=iacc(ki,:)
%     puu=[puu pMc+((ipu-1)*2*lKA+ipu-1)*nk1max];
%    end
%    sim=2*nk1max*length(iacc(ki,:));
%    xx=reshape(x(puu),sim,sim);
%    xm=xx^y;
%    P(puu)=xm;
%   end
%   P=P(Pust,Pust);
   P=IdeOo;
   sia=length(iacc);
   for ki=sia
     puu=PUrid{ki};
     sim=sqrt(length(puu));
     xx=reshape(x(puu),sim,sim);
     xm=xx^p;
     P(puu)=xm;
   end
  end
 elseif ifx==1

  P=IdeOo;
  for ki=1:lKA
   puu=pMu0u1+ki+(ki-1)*lKA*2;
   xx=reshape(x(puu),2,2);
   xm=xx^p;
   P(puu)=xm;
  end

  P=P(Pust,Pust);


 elseif ifx>=2
  isal=0;

  P=IdeOo;
 if ive==1
   if isal==1
   xp=x;
   yp=x;
   save xv xp yp
   ' powermat ve ', pausak
   end
  for ki=1:nk1max
   puu=pMu0u+ki+(ki-1)*2*lKA;
   xx=reshape(x(puu),nures,nures);
   xm=xx^p;
   P(puu)=xm;
  end
 else
   if isal==1
   xp=x;
   yp=x;
   save xn xp yp
   ' powermat nuo ', pausak
   end
  for ki=1:nk1max
   puu=pMu0u{ki};
   s=length(puu);
   nur=sqrt(s);
   xx=reshape(x(puu),nur,nur);
   xm=xx^p;
   P(puu)=xm;
  end
 end
  P=P(Pust,Pust);

 end

%toc,
%ifx
%pausak
