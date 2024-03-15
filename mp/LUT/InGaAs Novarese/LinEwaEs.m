clear all
close all

inew=1;  %2 con fattore BGR


%load Anders1060-6nmQw27In-Bar10nm_5Set19_

TV=[300];
%TV=[301];

for T=TV

if inew==1
 Ila=2;
 eval(['load Anders1060-6nmQw318In-Bar10nm_17Dic19renF',num2str(T)]')
elseif inew==0
 Ila=4;
 eval(['load Anders1060-6nmQw27In-Bar10nm_5Set19_',num2str(T)]')
elseif inew==2
 Ila=2;
 eval(['load Anders1060-6nmQw318In-Bar10nm_19Dic19ren',num2str(T)]') 
end

Densityv=N;
lambdavet=lambda;


hcel=6.62e-34*3e17/1.6e-19;    

hcl=ones(length(Densityv),1)*(hcel./(lambdavet));    

DE=repmat(EFspV,1,length(lambdavet));
DEL=(hcl-DE)/.025;


nsp=1./(1-exp(DEL));

fe=1e-12;

figure, plot(Densityv*fe,Es,'linewidth',2)
chold
plot(Densityv*fe,G.*nsp*0.8e10,'--')
ylim([0 8e13])
grid
xlabel(' N  (1e12/cm^2)')
ylabel('R_{sp} (1/s)')
pausak
figure, plot(Densityv*fe,Es./G,'linewidth',2)
xlabel(' N  (1e12/cm^2)')
ylabel('Effective n_{sp}')
ylabel('Effective n_{sp} x v_g')

pausak

Line



LaMat=ones(length(Densityv),1)*(.02*pi./(lambdavet*1e-9));    
gg=G./LaMat;  
nn=squeeze(Dep);

%subplot(144)
%plot(Densityv(1:end-1),(diff(nn)./diff(gg)).^1)



%ylabel('\alpha')
%xlabel('Carrier Dens')  

Line

end

