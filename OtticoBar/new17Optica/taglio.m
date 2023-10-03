%' passo da taglio ', keyboard
   if icut==1
    if iLP==1
     y=abs(ADouti);
    else
     lP=length(Pus);
%     pu1=[1:lP/2];
%     pu2=[lP/2+1:lP];
     pud1=[1:lP/2];
     pud2=[lP/2+1:lP];
     y=abs(ADouti(pud1))+abs(ADouti(pud2));
    end

%   if icut==1
    yp=diff(y);
    lyp=length(yp);
    xd=[1:lyp];
    xd0=[1:lyp+1];
    ypz=yp(1:lyp-1).*yp(2:lyp);
    fiypz=find(ypz<0)+1;
    ys=[diff(yp); 0];
    fima=find(ys(fiypz)<0);
    fimi=find(ys(fiypz)>0);
    yMa=y(fiypz);
    [YMA,fiMa]=max(yMa);
    fiMi=fiMa+1;
    if length(fiypz)<fiMi
     ipuma=length(KK);
    else
     ipuma=fiypz(fiMi);
    end
    if ifp>0
     figure, plot(xd0,y,'r',xd0(fiypz(fima)),y(fiypz(fima)),'wo',...
     xd0(fiypz(fimi)),y(fiypz(fimi)),'go',xd0(ipuma),y(ipuma),'c*')
     pausak
    end

   else
    ipuma=length(KK);
   end


 if iLP==1
  puAc=[1:ipuma];
  puAc2=[1+ipuma:length(KK)];
 else
  pu=[1:ipuma];
  puAc=[pu pu+length(KK)];
  pu=[1+ipuma:length(KK)];
  puAc2=[pu pu+length(KK)];
 end

