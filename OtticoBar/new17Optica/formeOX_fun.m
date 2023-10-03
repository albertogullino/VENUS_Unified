function [rv,fasid,fasud]=formeOX(R0x,DeFor,DeForC,NAz,dF)
  
  fi=linspace(0,2*pi,100001);
  cce=0;
  R00=1;
  %DeFor=.15;
  %DeForC=+.02;
  %NAz=3;
  %dF=pi/2;
  %dF=0;
  
  rocf=cce+R00*(1+DeFor*sin(NAz*fi)+DeForC*sin(2*NAz*fi)).*exp(j*fi);
  roc=R00*(1+DeFor*cos(NAz*(fi-dF))+DeForC*cos(2*NAz*(fi-dF)));
  rocf=cce+roc.*exp(j*fi);
  droc=-DeFor/NAz*sin(NAz*(fi-dF))-DeForC/(2*NAz)*sin(2*NAz*(fi-dF));

  ifi=0;
  if ifi==1
  figure, 
  plot(real(rocf),imag(rocf),'r','linewidth',2)
  axis equal
  grid
  pausak
  
  figure, plot(fi,roc,fi,droc+R00), pausak
  end
  %return
  
  fdcs=droc(1:end-1).*droc(2:end);
  fzder=find(fdcs<0 & droc(1:end-1)<0)+1;
  fzder1=fzder(1);
  Fio=fi([fzder1:end]);
  Fio=[Fio fi(1:fzder1-1)+2*pi];
  Fa=roc([fzder1:end 1:fzder1-1]);
  Fad=droc([fzder1:end 1:fzder1-1]);
  
  %Fa=roc;
  Fmi=min(Fa);
  Fma=max(Fa);
  rVet=linspace(Fmi,Fma,501);
  kinc=0;
  for krv=2:length(rVet)-1;
   kinc=kinc+1;
   r0=rVet(krv);
   fze=Fa-r0;
   fcs=fze(1:end-1).*fze(2:end);
   fizdu=find(fcs<=0);
   fizv=find(diff([fizdu fizdu(end)+10])>1);
   fiz=fizdu(fizv);
   deR=Fad(fiz);
   fizp=find(deR>=0);
   fizn=find(deR<0);
   if fizn(1)>fizp(1)
    FZi(kinc,1:length(fizp))=Fio(fiz(fizp));
    FZf(kinc,1:length(fizp))=Fio(fiz(fizn));
   else
    FZf(kinc,1:length(fizp))=Fio(fiz(fizp));
    FZi(kinc,1:length(fizp))=Fio(fiz(fizn));   
   end
  end
  rVp=rVet(2:end-1);

 if ifi==1
  figure, 
  plot(Fio,Fa,Fio,Fad+R00), pausak  
  
    figure, 
    plot(rVp,FZi,'g.',rVp,FZf,'r.'), pausak
  end

  fi=find(FZi~=0);
  FZid=FZi(fi);
  fi=find(FZf~=0);
  FZfd=FZf(fi);
  mm=min(min(FZid));
  MM=max(max(FZfd));
  dm=MM-mm;
  drpi=2*pi-dm;
  mm=mm-drpi/2;
  MM=MM+drpi/2;
  
  NP=301;
  FZia=linspace(0,Fmi,NP)';
  FZmN=ones(NP,size(FZi,2))*NaN;
  FZm=FZmN;
  FZm(:,1)=ones(size(FZia))*mm;
  FZp=FZmN;
  FZp(:,end)=ones(size(FZia))*MM;
  
  fasid=[FZm; FZi];
  fasud=[FZp; FZf];
  rv=R0x*[FZia; rVp'];
 if ifi==1
  figure, plot(rv,fasid,'g.',rv,fasud,'r.')
 end    
