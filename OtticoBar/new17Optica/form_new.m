icolo=0;
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
    lfi=100;
    fi=linspace(0,2*pi,lfi);
    ro=ones(1,lfi)*avero;
    romi=ro-dstr;
    roma=ro+dstr;
    rot=[roma(1) romi roma(1:length(roma))];
    fito=[0 fi -fi(1:length(roma))];
    roc1=rot.*exp(j*fito);
    roc2=roc1;
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
    romi=ro-dstr;
    roma=ro+dstr;
    rot=[roma(1) romi roma(1:length(roma))];
    fito=[0 fi -fi(1:length(roma))];
    roc1=rot.*exp(j*fito);
    roc2=roc1;
   elseif sha==3
    x0=linspace(bvero,avero+1e-6,25);
    ro=[x0 fliplr(x0) x0 fliplr(x0)];
    fi0=atan(sqrt(((x0/bvero).^2-1)./(1-(x0/avero).^2)));
    fi=[fi0 pi-fliplr(fi0) fi0+pi 2*pi-fliplr(fi0)];
    romi=ro-dstr;
    roma=ro+dstr;
%    rot=[roma(1) romi roma(1:length(roma)-1)];
%    fito=[0 fi -fi(1:length(roma)-1)];
    rot=[roma(1) romi roma(1:length(roma))];
    fito=[0 fi -fi(1:length(roma))];
    roc1=rot.*exp(j*fito);
    roc2=roc1;
   elseif (sha==4 | sha==5)
     npf=100;
     dmi=avero-dstr;
     dma=avero+dstr;
     fasi1=linspace(0,2*pi,npf);
     fasi2=linspace(0,2*pi,npf);
     mo1=linspace(1,1,npf);
     mo2=linspace(1,1,npf);
     mo=[dma mo1*dmi mo2*dma];
     fas=[0 fasi1 fasi2];
     if sha==4

          nv_sh=4;
          xt=fas;
          if istrumix==0
           rapax=bvero/avero;
          else
           rapax=b/a;
          end
          dae=(rapax-1)/2;
          ru=(1+dap*cos(nv_sh*xt)+dae*cos(2*xt));
          roc2=mo.*ru.*exp(j*fas);
     else
        Rvt=[];
        for nh=1:length(cce)
         roplo=cce(nh)+mo.*exp(j*fas);
%         Rvt=[Rvt NaN roplo];
         Rvt=[Rvt  roplo'];
        end
        roc2=Rvt.';

     end

   elseif sha==6

     X=pdi;
     D_stri=aloc/kcav0;
     d_stri=bloc/kcav0;
     if D_stri<0
      D_stri=-D_stri;
      orgra=pi/2;
     else
      orgra=0;
     end

     ns=round(2*X/D_stri);
     yv0=[-d_stri*ones(1,100) d_stri*ones(1,100)]/2;

     xcu=[];
     ycu=[];

     for k=1:ns
      Ym=(k-1-fix(ns/2))*D_stri;
      yvd=yv0+Ym;
      Xm=sqrt(X^2-(Ym)^2);
      xpar=linspace(-Xm,Xm,100);
      xvd=[xpar fliplr(xpar)];
      if orgra==0
       yv=yvd;
       xv=xvd;
      else
       xv=yvd;
       yv=xvd;
      end
      xcu=[xcu  [xv xv(1)]'];
      ycu=[ycu  [yv yv(1)]'];
     end

   end

   if length(roc2)>0
     xcu=real(roc2)';
     ycu=imag(roc2)';
     roc2=[];
   end
   zcu=ones(size(xcu));

%figure,
%fill3(xcu,ycu,zcu,'w')
%view(2), axis equal
%disp('form_new')
%keyboard
