% Sources normalized to their corresponding Thermal power

% Joule
figure
set(gcf,'pos',[631         454        1276         420])
subplot(121)
% surf(x,y,real(log10(Joule))),title('Joule (log)')
surf(x,y,real(log10(Joule/PJoule(iCurr-1)))),title('Joule (log)')
shading interp
view(2)
colorbar
xlabel('\rho, \mum')
ylabel('z, \mum')
set(gca,'FontSize',16,'FontName','Times new roman')
axis([0 8 114 118.5])
subplot(122)
% surf(x,y,Joule),title('Joule (lin)')
surf(x,y,Joule/PJoule(iCurr-1)),title('Joule (lin)')
shading interp
view(2)
colorbar
xlabel('\rho, \mum')
ylabel('z, \mum')
set(gca,'FontSize',16,'FontName','Times new roman')
axis([0 8 114 118.5])

% Ccap
figure
set(gcf,'pos',[631         454        1276         420])
subplot(121)
% surf(x,y,real(log10(Ccap))),title('Ccap (log)')
surf(x,y,real(log10(Ccap/Pcap(iCurr-1)))),title('Ccap (log)')
shading interp
view(2)
colorbar
xlabel('\rho, \mum')
ylabel('z, \mum')
set(gca,'FontSize',16,'FontName','Times new roman')
axis([0 8 114 118.5])
subplot(122)
% surf(x,y,Ccap,title('Ccap (lin)')
surf(x,y,Ccap/Pcap(iCurr-1)),title('Ccap (lin)')
shading interp
view(2)
colorbar
xlabel('\rho, \mum')
ylabel('z, \mum')
set(gca,'FontSize',16,'FontName','Times new roman')
axis([0 8 114 118.5])

% GR
figure
set(gcf,'pos',[631         454        1276         420])
subplot(121)
% surf(x,y,real(log10(GR))),title('GR (log)')
surf(x,y,real(log10(GR/Pnr(iCurr-1)))),title('GR (log)')
shading interp
view(2)
colorbar
xlabel('\rho, \mum')
ylabel('z, \mum')
set(gca,'FontSize',16,'FontName','Times new roman')
axis([0 8 114 118.5])
subplot(122)
% surf(x,y,GR),title('GR (lin)')
surf(x,y,GR/Pnr(iCurr-1)),title('GR (lin)')
shading interp
view(2)
colorbar
xlabel('\rho, \mum')
ylabel('z, \mum')
set(gca,'FontSize',16,'FontName','Times new roman')
axis([0 8 114 118.5])

% FCA
figure
set(gcf,'pos',[631         454        1276         420])
subplot(121)
% surf(x,y,real(log10(FCA))),title('FCA (log)')
surf(x,y,real(log10(FCA/Pfca(iCurr)))),title('FCA (log)')
shading interp
view(2)
colorbar
xlabel('\rho, \mum')
ylabel('z, \mum')
set(gca,'FontSize',16,'FontName','Times new roman')
axis([0 8 114 118.5])
subplot(122)
% surf(x,y,FCA),title('FCA (lin)')
surf(x,y,FCA/Pfca(iCurr)),title('FCA (lin)')
shading interp
view(2)
colorbar
xlabel('\rho, \mum')
ylabel('z, \mum')
set(gca,'FontSize',16,'FontName','Times new roman')
axis([0 8 114 118.5])