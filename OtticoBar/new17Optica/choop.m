load mat
D1=diag(Kplot1);
D2=diag(Kplot2);
Dm=D1+D2;
figure, plot(D1), hold on, plot(D2,'r')

lk=length(KK);
pu=1:lk;
pre=.05;
figure
Pus=[];
for nm=1:numodi
 pui=pu+(nm-1)*lk;
 Du=abs(Dm(pui));
 Du=Du/max(Du);
 der=diff([0; Du]);
 fiu=find(Du-pre>0 & der>0);
 pua=fiu(1):lk;
% plot(pu,Du,pu(pua),Du(pua),'r.')
% pausak
 Pus=[Pus pui(pua)];
end


 if ifp>=1
  lD=length(Dm);
  put=1:lD;
  figure,
  plot(put,Dm,put(Pus),Dm(Pus),'r.')
  [length(Dm) length(Pus)]
  map(Kplot2)
  map(Kplot2(Pus,Pus))
 end
