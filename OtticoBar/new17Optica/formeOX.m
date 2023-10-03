clear all
close all
  fi=linspace(0,2*pi,10001);
  
  iforma=input(' iforma = [0/1] ');
  if length(iforma)==0
   iforma=1;
  end
  
  cce=0;
  R00=1;
  vDeForP=[5 10 20 30];
  vDeForC=[-.02 .02];
  vDeForF=[0 20 40];
  %DeForC=+.0;
  %DeForS=-.001;
  %DeForS=.00;
  NAz=5;
  
  dF=pi/2;
  dF=0;
  
  figure
  for DeForP=vDeForP
  DeFor=DeForP/100;
  
  for DeForF=vDeForF
  
  DeForC=DeForF*DeFor/100;
  rocf=cce+R00*(1+DeFor*sin(NAz*fi)+DeForC*sin(2*NAz*fi)).*exp(j*fi);
  roc=R00*(1+DeFor*cos(NAz*(fi-dF))+DeForC*cos(2*NAz*(fi-dF)));
  rocf=cce+roc.*exp(j*fi);
  droc=-DeFor/NAz*sin(NAz*(fi-dF))-DeForC/(2*NAz)*sin(2*NAz*(fi-dF));

  
  plot(real(rocf),imag(rocf),'r','linewidth',2)
  axis equal
  a=axis;
  axis(a*1.2)
  grid
  title([' Par 1=',num2str(DeFor),' Par 2=',num2str(DeForC)])
  pausak
  
  end
  
  end
  
  if iforma==1
  return
  end
  figure, plot(fi,roc,fi,droc+R00), pausak

  %figure, plot(fi,roc,'r','linewidth',3), hold on
  %for mu=2:2:8
  % plot(fi,cos(mu*fi),fi,sin(mu*fi)), pausak
  %end

  
  %return
  
  fdcs=droc(1:end-1).*droc(2:end);
  fzder=find(fdcs<0)+1;
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
   fiz=find(fcs<0);
   deR=Fad(fiz);
   fizp=find(deR>=0);
   fizn=find(deR<0);
   if fizn(1)>fizp(1)
    FZi(kinc,:)=Fio(fiz(fizp));
    FZf(kinc,:)=Fio(fiz(fizn));
   else
    FZf(kinc,:)=Fio(fiz(fizp));
    FZi(kinc,:)=Fio(fiz(fizn));   
   end
  end
  rVp=rVet(2:end-1);

  figure, 
  plot(Fio,Fa,Fio,Fad+R00), pausak  
  
  figure, 
  plot(rVp,FZi,'g.',rVp,FZf,'r.'), pausak
  
  mm=min(min(FZi));
  MM=max(max(FZf));
  dm=MM-mm;
  drpi=2*pi-dm;
  mm=mm-drpi/2;
  MM=MM+drpi/2;
  
  NP=301;
  FZia=linspace(0,Fmi,NP)';
  FZmN=ones(NP,NAz)*NaN;
  FZm=FZmN;
  FZm(:,1)=ones(size(FZia))*mm;
  FZp=FZmN;
  FZp(:,end)=ones(size(FZia))*MM;
  
  fasid=[FZm; FZi];
  fasud=[FZp; FZf];
  rv=[FZia; rVp'];
  figure, plot(rv,fasid,'g.',rv,fasud,'r.')