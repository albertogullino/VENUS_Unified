asave=a;
if exist('b')
 bsave=b;
end
icolo=0;
fia=find(ar.a>0);
sha=shav.a(fia);
avero=ar.a(fia);

if sha>1 & sha<=3
 bvero=aral.y.a(fia);
else
 bvero=aral.y.a(fia);
 dap=aral.p.a(fia);
end
a=avero;
if exist('bvero')
 b=bvero;
end

if ~exist('ifalso')
 ifalso=0;
end

    fi=fian0;
    ro=ones(1,length(fi))*avero;

   if sha==1
    roc2=ro.*exp(j*fi);
   elseif sha==2
    fa=atan(a/b);
    f1=fi(find(fi<=fa));
    f2=fi(find(fi>fa & fi<=pi-fa));
    f3=pi-fliplr(f1);
    fu=[f1 f2 f3];
    fiu=[fu fu(2:end)+pi];
    r1=bvero*sqrt(1+tan((f1)).^2);
    r2=avero*sqrt(1+tan(pi/2-(f2)).^2);
    r3=bvero*sqrt(1+tan((f3)).^2);
    ru=[r1 r2 r3];
    ro=[ru ru(2:end)];
    roc2=ro.*exp(j*fiu);
   elseif sha==3
    x0=linspace(bvero,avero+1e-6,fix(length(fian0)/4));
    ro=[x0 fliplr(x0) x0 fliplr(x0)];
    fi0=atan(sqrt(((x0/bvero).^2-1)./(1-(x0/avero).^2)));
    fi=[fi0 pi-fliplr(fi0) fi0+pi 2*pi-fliplr(fi0)];
    roc2=ro.*exp(j*fi);
   elseif sha>=4
     mo=ro;
     if sha==4
          nv_sh=4;
          xt=fi;
          rapax=b/a;
          dae=(rapax-1)/2;
          ru=(1+dap*cos(nv_sh*xt)+dae*cos(2*xt));
          roc2=mo.*ru.*exp(j*fi);
     else
      disp(' in Borders: forma non prevista per il confinamento'), keyboard
        Rvt=[];
        for nh=1:length(cce)
         roplo=cce(nh)+mo.*exp(j*fas);
         Rvt=[Rvt  roplo'];
        end
        roc2=Rvt;

     end


   end


     xcu=real(roc2);
     ycu=imag(roc2);
     robor=abs(roc2);

a=asave;
if exist('bsave')
 b=bsave;
end
