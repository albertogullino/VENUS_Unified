%save fisal
%'salvato fis', keyboard

sp=6;
fsr=3;
if ~exist('STZ')
 STZ=1;
end
if ~exist('num_long')
 num_long=0;
end
if exist('iord_long')
 num_long=length(iord_long);
end
%sp=num_long*fsr;
%if num_long<3
% sp=3*fsr;
%end
lai=[0:fsr:sp]-sp/2;
lai=0;

lace=lambda0*1e9;
if exist('lambda_cen')
 lace=lambda_cen*1000;
end 
%lace=1650;
lave=lace+lai;
%'lave', keyboard
'lave', keyboard
clear Fivep Fivev gvetp gvetv lvetv lvetp
ifp_salvo=ifp;
for kla=1:length(lave)
 lambda0=lave(kla)*1e-9
  [Ksip,lambda,Fi,uLongp,uLong0p,Fap]=eiglmio(lambda0,uFunc,uF0,rel_pa,STZ);
  Gamp=uLongp;
  gpla = 2e-2*pi*rr/lambda*imag(Ksip)/Gamp;
 gvetp(kla)=gpla;
 lvetp(kla)=lambda*1e6;
 Fivep(:,kla)=Fi/max(Fi);
  [Ksip,lambda,Fi,uLongp,uLong0p,Fap]=eiglmio(lambda0,uFunc,uF0,rel_ve,STZ);
   Gamp=uLongp;
   gpla = 2e-2*pi*rr/lambda*imag(Ksip)/Gamp;
  gvetv(kla)=gpla;
  lvetv(kla)=lambda*1e6;
 Fivev(1:length(Fi),kla)=Fi/max(Fi);
 
end
rperm=rel_pa;

ifp=ifp_salvo;


[du,fiv]=min(gvetv);
lamv=lvetv(fiv)*1e-6;
[du,fip]=min(gvetp);
lamp=lvetp(fip)*1e-6;

%return

for kdir=[1 2]

if kdir==1
 gvet=gvetv;
 lvet=lvetv;
 Five=Fivev;
else
 gvet=gvetp;
 lvet=lvetp;
 Five=Fivep;
end

   [gs,iso]=sort(gvet);
   [las,isol]=sort(lvet*1000);
   dilas=diff([0 las]);
   fival=find(dilas>.1);
   dilasv=dilas(fival);
   FSR=mean(dilasv(2:end))
   laval=las(fival);
   gsoval=gvet(isol(fival));
   Fval=Five(:,isol(fival));

pulo=[1:length(gsoval)];   
%if iord_long>0   
 gsoval0=gsoval;
 laval0=laval;
   [gord,imi]=sort(gsoval);
   gsoval=gsoval(imi(pulo));
   laval=laval(imi(pulo));
   Fidis=Fval(:,imi(pulo));
   [du,ficond]=sort(laval);
   ficon=pulo;

      gso1=gsoval(ficon);
      laref=laval(ficon);
      lambdadu=laval(ficon)/1000;
      if kdir==1
       lav=lambdadu*1e-6;
       gav=gso1;
      else
       lap=lambdadu*1e-6;
       gap=gso1;
      end
x=hz;
if ifp==-10   
   figure, semilogy(laval0,gsoval0,'g.-',laval(ficon),gsoval(ficon),'ro'), pausak
      figure, plot(x,3*Fidis(:,ficon)),
      hold on, pco=plot(x,sqrt(real(rperm)),'b','linewidth',2);
' cont long', keyboard
      
end

%'fermo ', keyboard   

end   

miv=min(min(gav));
mip=min(min(gap));
ipol=1;   %parallela (p)
if mip>miv
 ipol=2;
end 

if ipol==1
 [mip,icol]=min(min((gap)));
 gap=gap(:,icol);
 lap=lap(:,icol);
 [gs,is]=sort(gav,2);  %ordino per righe
 ls=lav(is);
 gav=gs(:,1);
 lav=ls(:,1);
else
 [mip,icol]=min(min((gav)));
 gav=gav(:,icol);
 lav=lav(:,icol);
 [gs,is]=sort(gap,2);  %ordino per righe
 ls=lap(is);
 gap=gs(:,1);
 lap=ls(:,1);
end

if ifp==-10

%figure, semilogy(lave,gvetp,'ro',lave,gvetv,'go'), pausak
%figure, plot(lave,lvetp,'ro',lave,lvetv,'go'), pausak
figure, plot(lvetp,gvetp,'ro',lvetv,gvetv,'go'), 
hold on, plot(lap*1e6,gap,'rx','linewidth',5), 
hold on, plot(lav*1e6,gav,'gx','linewidth',5), 

pausak



%   figure, plot(x,rperm,'w',x,Five*3), 
%       title([' lambda_{res} = ',num2str(lambda),'  Gth = ',num2str(gpla)]), pausak
end      
      
