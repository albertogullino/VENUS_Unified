%========================================================================76
%clear
clear global
close all
colordef white
dbstop if error

addpath('Termico')
addpath('Ottico')
addpath('Dati')
addpath('out')

inewpat=1;
if inewpat==1
 rmpath('generageom')
 rmpath('Ottico')
 rmpath('Ottico\new17Optica')
 addpath('generageomBar')
 addpath('OtticoBar') 
 addpath('OtticoBar\new17Optica') 
else 
 rmpath('generageomBar')
 rmpath('OtticoBar\new17Optica') 
 addpath('generageom')
 addpath('Ottico\new17Optica')
end 



ilo=input(' Carica du1 ');
if length(ilo)==1

        load du1
        
end        

iCOMP=0;
        
mode.iTfig=0;
IPLOT=0;
        dTQW=T(mesh.inMQW{2})-T0;
        leQW=length(modePlot.nQW{end}{2});
        x=modePlot.x;
        xQW=modePlot.x(1:leQW);
        
        dTm=dTQW(1);
        FR=[1 2/3 1/2 1/3];
        clear T2 MT
        
        for kt=1:length(FR)
        [du,fi2]=min(abs(dTQW-dTm*FR(kt)));
        x2(kt)=xQW(fi2);
        T2(kt)=dTQW(fi2);
        end
        
        xQf=linspace(0,max(xQW),2001);
        
        figure, plot(xQW,dTQW,x2,T2,'ro'), pausak

Fattore=linspace(.1,1,11);
        
         for it=1:length(Fattore)
           mesh.fCondTer=Fattore(it);
          for iz=1:length(Fattore)
            mesh.fCondTerZ=Fattore(iz);
            [DeltaT,Tprec,PTherm,T_Contributi]=f_ThermicFun(Tprec,mesh,mode,StrTT,IPLOT);
%            'dopo f_Te', keyboard
                dTQW=DeltaT(mesh.inMQW{2});
                dTQWf=spline(xQW,dTQW,xQf);
                dTm=dTQW(1);                
%                Tprec=DeltaT+T0;
	        for kt=1:length(FR)
 	         [du,fi2]=min(abs(dTQWf-dTm*FR(kt)));
	         x2(kt)=xQf(fi2);
	         T2(kt)=dTQWf(fi2);
	        end
        
                figure(69),
                hold on
                plot(xQf,dTQWf,x2,T2,'ro'), 
                %pausak            
                
                Mx(it,iz,:)=x2;
                Mt(it,iz,:)=T2;
                MT(it,iz,:)=dTQWf;
                
          end
         end 
        
        Tm0=60;
        figure, semilogy(Fattore,Mt(:,:,1)),
        hold on, semilogy(Fattore,ones(size(Fattore))*Tm0,'r','linewidth',2),
        pausak
        figure, plot(Fattore,Mt(:,:,1))
            hold on, plot(Fattore,ones(size(Fattore))*Tm0,'r','linewidth',2),
            ylim(Tm0*[.2 1.5])

Fatf=linspace(Fattore(1),Fattore(end),101);

for kf=1:length(Fattore)
  Mf=Mt(:,kf,1);
  Tmf=spline(Fattore,Mf,Fatf); 
   [du,fi2]=min(abs(Tmf-Tm0));
   FaTr(kf)=Fatf(fi2);
  Mf=Mx(:,kf,2);
  Xd(kf)=spline(Fattore,Mf,Fatf(fi2));    
end

%fCondTer=.75;   % transverse thermal conducibility
%fCondTerZ=.6;   % Longitudinal thermal conducibility

FaZ0=0.6;
FaT0=0.75;
figure, plot(Fattore,FaTr,FaZ0,FaT0,'ro')
xlabel(' fCondTer Z')
ylabel(' fCondTer T')
            
            
            
        keyboard
