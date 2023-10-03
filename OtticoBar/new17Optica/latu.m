load punta

   Base=PAlens.Base;
   Alt=PAlens.Alt;
   Sth=PAlens.Sth;
   
   Ain=ailu(1);
   Afi=Base;
   if Ain>Afi
    'errore punta: rastremazione errata'
    keyboard
   end
   dBase=Afi-Ain;
   Avedis=[Ain:Sth:Afi];
   if Avedis(end)<Afi
    Avedis=[Avedis Afi];
   end
   lA=length(Avedis);
   Avedis=linspace(Ain,Afi,lA);
   dAve=diff(Avedis);
   Avedis=Avedis(2:end);
   nlA=length(dAve);
   Ztip=dAve/dBase*Alt;

adu=[ailu; ailu(end)];
dadu=diff(adu);
fiad=find(dadu~=0);
iin=1
clear dn an nn
for kf=1:length(fiad);
 uu=fiad(kf);
 puu=iin:uu
 iin=uu+1;
 dn(kf)=sum(dilu(puu)); 
 an(kf)=ailu(uu); 
 nn(kf,:)=nilu(uu,:); 
% kf, pausak
end

'fiad', keyboard
fiad0=find(dadu==0);
   
   ailusa=an';
   dilusa=dn';
   dilua=dn';
   dilua(1)=dn(1)+Ztip(1);
   ailu=[fliplr(Avedis)'; ailusa];
   dilu=[fliplr(Ztip)'; dilua];

figure, plot(ailu,cumsum(dilu),'.')   
   
   
   