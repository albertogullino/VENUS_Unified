%for kv=2:length(modePlot.ii_dd)         
% n2D=modePlot.Elqw(kv,:);
% p2D=modePlot.Hoqw(kv,:);
% TQW=modePlot.Tqw(kv,:);         
% np2D=n2D.*p2D-n2Di.*p2Di;
% tauleak=tau_lea*f_RaisedCosine(TQW,Tstart,DeT,taubottom,tautop);              
% denLeakage2D = tauleak.*(n2D+n2Di)+tauleak.*(p2D+p2Di);
% Rlea = NumQW*qel*np2D./denLeakage2D;
% Nlea(kv)=2*pi*trapz(x,x.*Rlea)/Ccap(kv);
% 'ver;,', keyboard
%end 

% I=mode.ii_dd*1000;
% V=mode.vv_dd;
% P=sum(mode.Pst_dd,1)+mode.Psp_dd;


Rec=mode.IntRec;
Ccap=-sum(mode.IntCcapN,2);
Rnor=diag(1./Ccap)*Rec;

Rnor(:,5)=Rnor(:,3)*(1-mode.Fat_regeneration);

he=figure;
set(he,'pos',[349         398        1525         525 ])
subplot(121)
         semilogy(mode.ii_dd*1000,abs(Rnor),'LineWidth',1.5), 
%         legend('IntSRH','IntRad','IntAug','IntStim','IntLea','location','best')
         legend('IntSRH','IntRad','IntAug','IntStim','location','best')
         ylabel('Recombinations normalized to Injection')
         xlabel('Current (mA)')
         ylim([1e-3 1])
         
subplot(122) 
%Etai=Rnor(:,4)./sum(Rnor(:,[1 2 3 5]),2);
Etai=Rnor(:,4)./sum(Rnor,2);
Pott=sum(mode.Pst_dd,1);
corr=mode.ii_dd;
volt=mode.vv_dd;
Watt=corr.*volt*1000;
Etawp=Pott./Watt;
Eta=[Etai Etawp'];
         plot(mode.ii_dd*1000,Eta,'LineWidth',1.5), 
         legend('Internal','Wallplug','location','best')
         title('Efficiency')
         xlabel('Current (mA)')
         
         [du,ipm]=max(Pott);
         [du,ipw]=max(Etawp);
         fci=figure;
          set(fci,'pos',[456         343        1294         525])
         subplot(121)
         plot(Watt,Eta,'LineWidth',1.5)
                  hold on, plot(Watt(ipm)*ones(1,2),Eta(ipm,:),'ko','LineWidth',2)
                  hold on, plot(Watt(ipw)*ones(1,2),Eta(ipw,:),'go','LineWidth',2)
	          legend('Internal','Wallplug','location','best')
	          ylabel('Efficiency, %')
	          grid
	          xlim([0 40])
         xlabel('Electrical Power, mW')
                  subplot(122)
	          plot(Watt,Pott,'LineWidth',1.5)
	                   hold on, plot(Watt(ipm),Pott(ipm),'ko','LineWidth',2)
	                   hold on, plot(Watt(ipw),Pott(ipw),'go','LineWidth',2)
	 	          ylabel('Light, mW')
         xlabel('Electrical Power, mW')
         	          xlim([0 40])
         grid