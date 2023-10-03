clear all
close all
load sa

iplo=0;
%Par=Ps.Par;
%mazim=Par.mazim;
%SogDen=Par.SogDen;

nubesu=4;
mazim=4;
SogDen=0;

an1=acos(SogDen);
an2=acos(-SogDen);

pea=2*pi/mazim;

anv0=[an1 an2]/mazim;
anv=[anv0 pea-(anv0)];
sogv=[SogDen -SogDen SogDen -SogDen];


an_vet=[];
so_vet=[];
for k=1:mazim
 an_vet=[an_vet anv+(k-1)*pea];
 so_vet=[so_vet sogv];
end

an_veti=an_vet(1:2:end-1);
an_vetu=an_vet(2:2:end);
[an_veti' an_vetu']

firad=linspace(0,2*pi,501);

if iplo==1
figure, plot(firad,cos(mazim*firad),an_vet,so_vet,'ro'), 
hold on, plot(firad,sin(mazim*firad),'r',an_vet,sin(mazim*an_vet),'go'),
plot(firad,ones(size(firad))*SogDen,'w--',firad,-ones(size(firad))*SogDen,'w--')
pausak

figure, plot(an_veti,sin(mazim*an_veti),'go',an_vetu,sin(mazim*an_vetu),'ro'), pausak
figure, plot(an_veti,cos(mazim*an_veti),'go',an_vetu,cos(mazim*an_vetu),'ro'), pausak
end

IcR=(mazim*diff(anv(1:2)) +1/(4*mazim)*sum(sin(mazim*an_vetu)-sin(mazim*an_veti)))/pi;
IsR=(mazim*diff(anv(1:2)) -1/(4*mazim)*sum(sin(mazim*an_vetu)-sin(mazim*an_veti)))/pi;

IR0=mazim*diff(anv(1:2));

ic=0;
pufi=2:24
clear Ic Is
for nufi=pufi
 ic=ic+1;
 Ic(ic)=1/(4*nufi)*sum(cos(nufi*an_vetu)-cos(nufi*an_veti));
 Is(ic)=-1/(4*nufi)*sum(cos(nufi*an_vetu)-cos(nufi*an_veti));
end

if iplo==1
figure, plot(pufi,Ic,pufi,Is,'ro'),  
title('Ic, Is')
pausak
figure, plot(pufi,Ic+IR0,pufi,Is+IR0,'ro'),  
title('IcTot, IsTot')
xlabel(' Ordine azimutale')
pausak
end

    lfi=31;
    fi=linspace(an_veti(1),an_vetu(1),lfi);
    ro=linspace(adis,bdis,11)/kcav;
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
 ftot=[ftot ftot1'+(kp-1)*pea/2];
 rtot=[rtot rtot1' ];
end 

    lfi=301;
    fi=linspace(0,2*pi,lfi);
if iplo==1

figure, plot(rtot.*exp(j*ftot))
hold on
plot(ro(1)*exp(j*fi),'g--')
plot(ro(end)*exp(j*fi),'r--'), axis equal
pausak

end
%figure, plot(firad,cos(mazim*firad),an_vet,so_vet,'ro'), pausak


shape='sha_vortex';

kmat_Ve2