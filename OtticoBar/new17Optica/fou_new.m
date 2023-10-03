%' entro fou ', pausak
%keyboard
%       if ~exist('sha_plo')
        sha_plo=sha;
%       end
global cce
sha_Cug=sha;
if DOE==1
 sha_Cug=8;
end

iset=0;
ipro=1;
icolo=0;
iarc(1,1)='w';
iarc(2,1)='y';
iarc(3,1)='r';
iarc(4,1)='g';
iarc(5,1)='c';
iarc(6,1)='b';
avero=aloc/kcav0;
if exist('bloc')
 bvero=bloc/kcav0;
end
if ~exist('ifalso')
 ifalso=0;
end
s_line=.2;  %micron
dstr=s_line/2;


   if sha==1
    lfi=101;
    fi=linspace(0,2*pi,lfi);
    ro=ones(1,lfi)*avero;
%    romi=ro-dstr;
%    roma=ro+dstr;
%    rot=[roma(1) romi roma(1:length(roma))];
%    fito=[0 fi -fi(1:length(roma))];
%    roc1=rot.*exp(j*fito);
    roc1=ro.*exp(j*fi);
    roc2=roc1;
   elseif sha==11
    if SPP==0
    lfi=31;
    fi=linspace(an_veti(2),an_vetu(2),lfi);
    ro=linspace(aloc,bloc,11)/kcav;
    f1=fi(1)*ones(size(ro));
    r1=ro;
    f2=fi;
    r2=ro(end)*ones(size(fi));
    f3=fi(end)*ones(size(ro));
    r3=fliplr(ro);
    f4=fliplr(fi);
    r4=ro(1)*ones(size(fi));
    ftot1=[f1 f2 f3 f4];
    rtot1=[r1 r2 r3 r4];
	ftot=[];
	rtot=[];
	for kp=1:mazim*2
	 ftot=[ftot ftot1'+kp*pea/2];
	 rtot=[rtot rtot1' ];
	end 
     roc1=rtot.*exp(j*ftot);
     roc2=roc1;	   
     iset=1
        sr=size(roc2);
%        roc2=[ones(sr(1),1)*sha_Cug ones(sr(1),1)*sr(2) roc2];
        roc2=[ones(1,sr(2))*sr(2); ones(1,sr(2))*length(rtot); roc2];
        xcu=real(roc2);
        ycu=imag(roc2);
        roc2=[];
    else  %SPP
      lfi=51;
      fi=linspace(an_veti,an_vetu,lfi)';
      
      roc1=[aloc/kcav.*exp(j*fi); bloc/kcav.*exp(j*flipud(fi)); aloc/kcav.*exp(j*fi(1))];
      sr=size(roc1);
      roc2=[ones(1,sr(2))*sr(2); ones(1,sr(2))*length(roc1); roc1];
      xcu=real(roc2);
      ycu=imag(roc2);
      iset=1;
      roc2=[];      
%     'founow', keyboard      
    end
   elseif sha==2
    fa=atan(a/b);
    f1=linspace(0,fa,10);
    f2=linspace(fa,pi-fa,30);
    f3=pi-fliplr(f1);
    fu=[f1 f2 f3];
    fi=[fu fu+pi];
    r1=bvero*sqrt(1+tan(f1).^2);
    r2=avero*sqrt(1+tan(pi/2-f2).^2);
    r3=bvero*sqrt(1+tan(f3).^2);
    ru=[r1 r2 r3];
    ro=[ru ru];
%    romi=ro-dstr;
%    roma=ro+dstr;
%    rot=[roma(1) romi roma(1:length(roma))];
%    fito=[0 fi -fi(1:length(roma))];
%    roc1=rot.*exp(j*fito);
%    ro=ru;
    roc1=ro.*exp(j*fi);
    roc2=roc1;
   elseif sha==3
    x0=linspace(bvero,avero+1e-6,25);
    ro=[x0 fliplr(x0) x0 fliplr(x0)];
    fi0=atan(sqrt(((x0/bvero).^2-1)./(1-(x0/avero).^2)));
    fi=[fi0 pi-fliplr(fi0) fi0+pi 2*pi-fliplr(fi0)];
    roc1=ro.*exp(j*fi);
    roc2=roc1;
   elseif (sha==4 | sha==5)
     npf=100;
     dmi=avero;
     fasi1=linspace(0,2*pi,npf);
     mo1=linspace(1,1,npf);
     mo=mo1*dmi;
     fas=fasi1;
     if sha==4
%     ' sha = 4 forma', keyboard
         nv_sh=4;
         xt=fas;
        if cce~=0
          Pshid{1}=sha;
          Pshid{2}(1)=bloc/kcav0;
          Pshid{3}(1)=aloc/kcav0;
          Pshid{4}(1)=dap;
          Pshid{5}(1)=cce;
          roc2=sha_gen(Pshid,xt,1);
          roc2=roc2.';
        else  
          if istrumix==0
           rapax=bvero/avero;
          else
           rapax=b/a;
          end
          dae=(rapax-1)/2;
          rapax=bvero/avero;
          dae=-(1-rapax)/(1+rapax)*(1+dap);
          R00=bvero/(1+dae+dap);

          ru=R00*(1+dap*cos(nv_sh*xt)+dae*cos(2*xt));
          roc2=mo/avero.*ru.*exp(j*fas);
       end %cce
       if isfield(Ps,'LatiPG')==1
         NAz=Ps.LatiPG;
         R00=aloc/kcav; 
         
         if abs(aiat-R00)<1e-14
          DeForp=Ps.aDeFor;
          DeForCp=Ps.aDeForC;
         else
          DeForp=Ps.DeFor;
          DeForCp=Ps.DeForC;
         end 
         dF=Ps.dF;
 
   
         DeFor=DeForp/R00;         
         DeForC=DeForCp/R00;         
         roc=R00*(1+DeFor*cos(NAz*(xt-dF))+DeForC*cos(2*NAz*(xt-dF)));
         roc2=cce+roc.*exp(j*xt);
        end 
%          'controlla foubu',          keyboard
     else
       % Rvt=[];
       % for nh=1:length(cce)
       %  roplo=cce(nh)+mo.*exp(j*fas);
       %  Rvt=[Rvt  roplo'];
       % end


        lfip=length(find(real(cce)>=0 & imag(cce)>=0));

        if lfip==length(cce)
         ccep=cce;
         Pshp=Pshi;
         icon=length(cce);
         lfip=(find((real(cce)>0 & imag(cce)==0) |(real(cce)==0 & imag(cce)>0) ));
         if length(lfip)>0
          ccep=[ccep; -cce(lfip)];
          for kc=lfip'
           icon=icon+1;
           for kP=1:5
            if kP<5
             Pshp{kP}(icon)=Pshi{kP}(kc);
            else
             Pshp{kP}(icon)=-Pshi{kP}(kc);
            end
           end
          end
         end

         lfip=(find(real(cce)>0 & imag(cce)>0));
         if length(lfip)>0
          cced=[-cce(lfip); conj(cce(lfip));  -conj(cce(lfip))];
          ccep=[ccep; cced];
          for kc=lfip'
           for kPi=1:3
            icon=1+icon;
            for kP=1:5
             if kP<5
              Pshp{kP}(icon)=Pshi{kP}(kc);
             else
              if kPi==1
               Pshp{kP}(icon)=-Pshi{kP}(kc);
              elseif kPi==2
               Pshp{kP}(icon)=conj(Pshi{kP}(kc));
              elseif kPi==3
               Pshp{kP}(icon)=-conj(Pshi{kP}(kc));
              end
             end
            end

           end
          end
         end

        else
         ccep=cce;
         Pshp=Pshi;
        end

        for nh=1:length(ccep)
         roplo=sha_gen(Pshp,fas,nh);
         Rvt(1,nh)=sha;
         Rvt(2,nh)=length(roplo);
         Rvt([1:length(roplo)]+2,nh)=roplo;
         Rvt(length(roplo)+3,nh)=ndis;
        end

        roc2=Rvt.';
        xcu=real(roc2)';
        ycu=imag(roc2)';
        iset=1;

     end

   elseif sha==6

ipro=0;
%' fou_new', keyboard

if ipro==0

X=P.Rx;
Y=P.Ry;
Shap=P.shape;
D_stri=P.D;
d_stri=P.d;
r_add=P.Radd;
if P.orien~0
 orgra=pi/2;
 Xlo=X;
 Ylo=Y;
else
 orgra=0;
 Xlo=Y;
 Ylo=X;
end

 ns=round(2*Xlo/D_stri);
if P.shif~0
 shiret=-D_stri/2;
else
 shiret=0;
end
 if shiret==0
  if is_even(ns)==1
   ns=ns+1;
  end
 else
  if is_even(ns)==0
   ns=ns+1;
  end
 end
' parametri '
'[ ns Xlo Ylo D_stri d_stri Shap]'
[ ns Xlo Ylo D_stri d_stri]
if d_stri<.1
 return
end
%' qui ', keyboard

   switch Shap
    case 'circle'
     Rma=Y;
     cas=0;
    case {'square','rectangle'}
     cas=1;
     Cmae=(ns-1)/2*D_stri+d_stri/2;

     if orgra==0
      Rma=sqrt(Cmae^2+X^2);
     else
      Rma=sqrt(Cmae^2+Y^2);
     end
   end

npr=101*fix(Rma/D_stri);
npr=500;
rv=linspace(0,Rma,npr);
%pausak
%keyboard


  yv0=d_stri/2*[-1 1];
  npf=50;
  rtot=[];
  rtoti=[];
  for k=1:ns
   Ym=(k-1-fix(ns/2))*D_stri;
   yvd=yv0+Ym;
   aas=yvd/Rma;
   fias=find(abs(aas)<=1);

   fiv=asin(aas(fias));
   if length(fiv)~=0
    if length(fiv)==1
     if fiv<0
      five=linspace(pi-fiv,2*pi+fiv,2*npf);
     else
      five=linspace(fiv,pi-fiv,2*npf);
     end
     five=five+orgra;
     rfiu=Rma*exp(j*five);
     fiex=Rma*cos(five([1 end]));
     rfiix=linspace(fiex(1),fiex(2),2*npf);
     rfiiy=ones(size(rfiix))*Rma*sin(five(1));
     rfii=rfiix+j*rfiiy;
     rtoti=[rfii rfiu];
     rtot=[rtot; rtoti];
    else
     five1=linspace(fiv(1),fiv(2),npf)+orgra;
     five2=linspace(pi-fiv(2),pi-fiv(1),npf)+orgra;
     rfiu1=Rma*exp(j*five1);
     rfiu2=Rma*exp(j*five2);
     fiex=Rma*cos([five1(end) five2(1)]);
     rfiix=linspace(fiex(1),fiex(2),npf);
     rfiiy=ones(size(rfiix))*Rma*sin(five1(end));
     rfii1=rfiix+j*rfiiy;

     fiex=Rma*cos([five2(end) five1(1)]);
     rfiix=linspace(fiex(1),fiex(2),npf);
     rfiiy=ones(size(rfiix))*Rma*sin(five1(1));
     rfii2=rfiix+j*rfiiy;
     rtoti=[rfiu1 rfii1 rfiu2 rfii2];
     rtot=[rtot; rtoti];
    end
   end
%    plot3(real(rtot.'),imag(rtot.'),ones(size(rtot.')),'w'), view(2), hold on

%   pausak
  end

%    xtot=real(rtot.');
%    ytot=imag(rtot.');
%    figure, plot(xtot,ytot)
%    keyboard
    roc2=rtot;
    if r_add>0
     fN=find(abs(rtot)<r_add);
     roc2(fN)=NaN;
%     'roc'
%     fiang=linspace(0,2*pi,length(roc2))';
%     roc2=[roc2  r_add*exp(j*fiang)];
     fiang=linspace(0,2*pi,length(roc2));
     roc2=[roc2;  r_add*exp(j*fiang)];

%     fiang=linspace(0,2*pi,length(roc2));
%     roc2=[roc2';  r_add*exp(j*fiang)];
%     keyboard
    end
% ' QUI ', keyboard
else %ipro
%     X=pdi;
%     D_stri=aloc/kcav0;
%     d_stri=bloc/kcav0;
%     if D_stri<0
%      D_stri=-D_stri;
%      orgra=pi/2;
%     else
%      orgra=0;
%     end
%     if d_stri<0
%      d_stri=-d_stri;
%      shiret=-D_stri/2;
%     else
%      shiret=0;
%     end

     X=P.Rx;
     Y=P.Ry;
     D_stri=P.D;
     d_stri=P.d;
     r_add=P.Radd;
     Shap=P.shape;
     Rma=Y;

     if P.orien~0
      orgra=pi/2;
      Xlo=X;
      Ylo=Y;
     else
      orgra=0;
      Xlo=Y;
      Ylo=X;
     end

     if P.shif~0
      shiret=-D_stri/2;
     else
      shiret=0;
     end

     Xma=max([X Y]);
     ns=round(2*Xlo/D_stri);
     if shiret==0
      if is_even(ns)==1
       ns=ns+1;
      end
     else
      if is_even(ns)==0
       ns=ns+1;
      end
     end

     lfi=50;
     yv0=[-d_stri*ones(1,lfi) d_stri*ones(1,lfi)]/2;

     xcu=[];
     ycu=[];

     for k=1:ns
      Ym=(k-1-fix(ns/2))*D_stri+shiret;
      yvd=yv0+Ym;
      switch Shap
       case 'circle'
        dif=Xlo^2-(Ym)^2;
       case {'square','rectangle'}
        dif=Ylo^2;
      end
      if dif>0
       Xm=sqrt(dif);
       xpar=linspace(-Xm,Xm,lfi);
       xvd=[xpar fliplr(xpar)];
       if orgra==0
        yv=yvd;
        xv=xvd;
       else
        xv=yvd;
        yv=xvd;
       end
       if r_add>0
        fiad=find(xv.^2+yv.^2<r_add^2);
        xv(fiad)=NaN;
        yv(fiad)=NaN;
       end
       xcu=[xcu  [sha_Cug length(xv)+1 xv xv(1)]'];
       ycu=[ycu  [sha_Cug length(xv)+1 yv yv(1)]'];
      end
     end
    if r_add>0
     fi=linspace(0,2*pi,length(xv)+1);
     cirv=r_add*exp(j*fi);
       xcu=[xcu  [sha_Cug length(fi) real(cirv)]'];
       ycu=[ycu  [sha_Cug length(fi) imag(cirv)]'];
    end
    iset=1;
   end
end %ipro


%   if length(roc2)>0
   if iset==0
     sr=size(roc2);
     roc2=[ones(sr(1),1)*sha_Cug ones(sr(1),1)*sr(2) roc2];
     xcu=real(roc2)';
     ycu=imag(roc2)';
     roc2=[];
   end

%'qui forma', keyboard

   zcu=ones(size(xcu));

   Prele=10;
   if exist('Rel')
    if Rel>0
     Prele=abs(aloc/kcav0-Rel);
    end
   end
   condin=(ndis==1 | ndis==length(adis) | length(aloc)<20);
   if DOE==1 | sha_plo==7 | sha==11
    condin=1;
   end
%   condin=1
   if (ifp==-10 | ifp>=0) & (condin | Prele<.01)
%     if exist('figsha')==0
      figsha=figure;
   %  else
    %  figure(figsha);
    % end
     xcui=xcu;
     ycui=ycu;
     zcui=zcu;
     fczu=1;
      si=size(xcui);
      for jj=1:si(2)
       npo=xcui(2,jj);
       pux=3:2+npo;
       xcud=xcui(pux,jj);
       ycud=ycui(pux,jj);
       zcud=zcui(pux,jj);
       if ~exist('sha_plo')
        sha_plo=sha;
       end
%       'shaplo', keyboard
       if sha_plo~=5
        co='w';
        if sha==11
         co='r';
        end
%        'qui', keyboard
        ph=plot3(xcud,ycud,zcud,co);
        hold on
        if DOE==0
        set(ph,'linewidth',3)
        else
        set(ph,'linewidth',1)
        end
%        'passo fow', keyboard
        if sha==11 & exist('ndise')

         if ndise==1
          set(ph,'color',[0,1,0])
         else
          set(ph,'linestyle','--')         
         end
        end 

        if sha_plo==6
         set(ph,'linewidth',1)
        elseif sha_plo==7
         set(ph,'linewidth',0.5,'linestyle','--')
        end
%        pausak
       else
        lab=xcui(3+npo,jj);
%        if lab==1
         ph=plot3(xcud,ycud,zcud*fczu,iarc(lab));
%        elseif lab==2
%         ph=plot3(xcud,ycud,zcud*fczu,'r');
%        elseif lab==3
%         ph=plot3(xcud,ycud,zcud*fczu,'y');
%        end
        set(ph,'linewidth',1)
        hold on
       end
        view(2), axis equal
      end
    disp('form_new')
%    pausak
    drawnow
      if ifp>=0
       grid
       pausak
      end
%      if ifp==-10 & numodiacc*length(KK)<200
%       pausak
%      end
     if sha==max(shailoop) & ifp~=-4
      grid on
%      keyboard
       if numodiacc*length(KK)<200 & ifr==1
        if ndis==1 | ndis==length(aloc)
         pausak
        end
       end
     end

   end
%sha

%figure,
%plot3(xcu(2:end,:),ycu(2:end,:),zcu(2:end,:),'w')
%view(2), axis equal
%disp('fou_new'), pausak
%keyboard
