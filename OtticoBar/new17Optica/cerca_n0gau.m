%iplo=1

clear Gaf Taf Gau

Dpe=0;
Dpm=0;

NPF=10;
irefv=20;
NPF0=50;
 nmi=min(real([n1 n2]));
 nma=max(real([n1 n2]))*1.3;
 nmai=max(abs(imag([n1 n2])));
 if nmai>1
  nma=5.5;
  nstart=.1;
 else
  nstart=nmi+.01;
  nmai=3;
 end 
 
% 'quio', keyboard
if icomp==1
 enner=linspace(nstart,nma,NPF0);
% enner=linspace(nmi,nma,NPF0);
% enner=linspace(nmi,nma,NPF0);
else
 nme=(n_ve+n_pa)/2;
 nmi=min([n1 n2]);
 nma=max([n1 n2]);
 ennem=linspace(nmi/2,nme,NPF0);
 ennee=linspace(nme,nma+1,NPF0);
%  ennee=linspace(.1,4,NPF0);
%   ennem=linspace(.1,4,NPF0);
 enner=ennee;
end


enner0=enner;
str=diff(enner(1:2));
strm=str;
if icomp==1
 ennei=linspace(-nmai,.1,40);
% ennei=-linspace(0,3,50);
 sti=diff(ennei(1:2));
else 
 ennei=0;
 sti=0;
end
 for kl1=1:length(enner)
 for kl2=1:length(ennei)
 nif=enner(kl1)+i*ennei(kl2);
  [gd,du,tp,du,gdu]=gaemt(0,0,lami,thick,nif,Lpm,npm,npair,r_inint,rr,0,lu,nu,r_in);
   Taf(kl1,kl2)=[tp];
   Gaf0(kl1,kl2)=[gd];
   Gau(kl1,kl2)=[gdu];
 end
 end
 
if icomp==1 
 if itetm==1 | itetm==3
  Gace=Gac(1);
  Gaf=Gaf0+Taf.^2.*Gace./(1-Gau*Gace);
  Dpe=Gaf-Glae;
  [du,i1v]=(min(abs(Dpe)));
  if du(1)<=0
   i2=1;
  else 
   [du,i2]=(min(du));
  end 
  ire=i1v(i2);
  iie=i2;
  nf1=enner(i1v(i2))+i*ennei(i2);
 if iplo==1
  figure, plot(enner,abs(Dpe),real(nf1),0,'wo')
  'E0c', pausak
  %keyboard
 end
  else
   nf1=0;
  end 
  if itetm==2 | itetm==3
   if length(Gac)==2
    Gacm=Gac(2);
   else 
    Gacm=Gac;
   end 
   Gaf=Gaf0+Taf.^2.*Gacm./(1-Gau*Gacm);
   Dpm=Gaf-Glam;
   [du,i1v]=(min(abs(Dpm)));
   [du,i2]=(min(du));
   irm=i1v(i2);
   iim=i2;
   nf2=enner(i1v(i2))+i*ennei(i2);
   Dpun0=Dpm(irm,iim);
  else
   nf2=0;
  end 
 if iplo==1
  figure, plot(enner,abs(Dpm),real(nf2),abs(Dpun0),'wo')
  'H0c', pausak
  %keyboard
 end
else  % icomp
 if itetm==1 | itetm==3
  Dpe=Gaf-Glae;
  [du,i1v]=(min(abs(Dpe)));
  [du,i2]=(min(du));
  ire=i1v(i2);
  iie=i2;
  nf1=enner(i1v(i2))+i*ennei(i2); 
  if iplo==1
  % figure, plot(enner,abs(Gaf-Glae),nf1,0,'wo')
   figure, plot(enner,abs(Dpe),nf1,0,'wo')
   'E0', keyboard
  end
 else
  nf1=0;
 end

 if itetm==2 | itetm==3
  enner=ennem;
  for kl1=1:length(enner)
  for kl2=1:length(ennei)
  nif=enner(kl1)+i*ennei(kl2);
  [gd,du,tp,du,gdu]=gaemt(0,0,lami,thick,nif,Lpm,npm,npair,r_out,rr,0,lu,nu,r_inint);
   if iga==0
    Ga_p= gd;
   else
     if length(Gac)==1
      Gacm=Gac;
     else
      Gacm=Gac(2);
     end  
    Ga_p=gd+tp^2*Gacm/(1-gdu*Gacm);
   end   
    Gaf(kl1,kl2)=[Ga_p];
  end
  end
  Dpm=Gaf-Glam;
  ennem0=ennem;
  [du,i1v]=(min(abs(Dpm)));
  [du,i2]=(min(du));
  irm=i1v(i2);
  iim=i2;
  nf2=enner(i1v(i2))+i*ennei(i2);
  if i1v==NPF0
   ennem=linspace(nme,nma,NPF0);
   strm=diff(ennem(1:2));
   enner=ennem;
   clear Gaf
   for kl1=1:length(enner)
    for kl2=1:length(ennei)
     nif=enner(kl1)+i*ennei(kl2);
     [gd,du,tp,du,gdu]=gaemt(0,0,lami,thick,nif,Lpm,npm,npair,r_out,rr,0,lu,nu,r_inint);
     if iga==0
      Ga_p= gd;
     else 
      Ga_p=gd+tp^2*Gac/(1-gdu*Gac);
     end   
     Gaf(kl1,kl2)=[Ga_p];
    end
   end
   Dpm=Gaf-Glam;
   ennem0=ennem;
   [du,i1v]=(min(abs(Dpm)));
   [du,i2]=(min(du));
   irm=i1v(i2);
   iim=i2;
   nf2=enner(i1v(i2))+i*ennei(i2);
  end
 
  if iplo==1
 %  figure, plot(enner,abs(Gaf-Glam),nf2,0,'wo')
   figure, plot(enner,abs(Dpm),nf2,0,'wo')
   'H0', keyboard
   end
  else
   nf2=0;
  end

end  % icomp
Dpe0=Dpe;
Dpm0=Dpm;


%' dppo', keyboard

stre=str;
stie=sti;
stim=sti;

%' sono quoi', keyboard
for lref=1:irefv
if itetm==1 | itetm==3
refegau
end
if itetm==2 | itetm==3
refmgau
end
end

if iplo==1
 figure, 
if itetm==1 | itetm==3
subplot(211)
 plot(enner0,abs(Dpe0),ennere(ire),abs(Dpee(ire,iie)),'wo',npBV,0,'ro'), 
 title([' Procedura determinazione  ref-index TE: n_{eff}=   ',num2str(nf1)])
end
if itetm==2 | itetm==3
  subplot(212)

 if icomp==0
  plot(ennem0,abs(Dpm0),ennerm(irm),abs(Dpem(irm,iim)),'wo',nvBV,0,'ro'), 
 else
 
  plot(enner0,abs(Dpm0),ennerm(irm),abs(Dpem(irm,iim)),'wo',nvBV,0,'ro'),  
 end
end
 title([' Procedura determinazione  ref-index TM: n_{eff}=   ',num2str(nf2)])

 pausak
end

return

if iplo==1
if itetm==1 | itetm==3
 figure, plot(enner0,abs(Dpe0),ennere(ire),abs(Dpee(ire,iie)),'wo'), 
 pausak
end
if itetm==2 | itetm==3
 if icomp==0
  figure, plot(ennem0,abs(Dpm0),ennerm(irm),abs(Dpem(irm,iim)),'wo'), 
 else
  figure, plot(enner0,abs(Dpm0),ennerm(irm),abs(Dpem(irm,iim)),'wo'),  
 end
 pausak
 end
end