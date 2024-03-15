
indLivello =1
xv = reshape(xvCtmp(:,indLivello),4,mesh.nn);
figure,plot(mesh.x,abs(sum(xv,1)))
hold on
plot(mesh.x,abs(xvCtmps(:,1)),'k')

% indLivello = 1
% 
% xv = reshape(xvVtmp(:,indLivello),4,mesh.nn);
% figure,plot(mesh.x,sum(xv,1))
% 



kxvet = linspace(0.000001,0.2,1001);
% kyvet = zeros(size(kxvet));

Bands = zeros(4,length(kxvet));

for indk = 1:length(kxvet)
    xmol = 0;%quellochevuoi
    strain = 0;%quell'altrochevuoi
        [H0,~,~,~]=assem_kp44_InGaAs(kxvet(indk),0,0,0);
       % [H0,~,~,~]=assem_kp44_AlGaAs(kxvet(indk),0,xmol);
    [V,E] = eig(H0);
    Bands(:,indk) = diag(E);
end
qel=1.6021766208e-019; % Elementary charge, C
Eg = 0.354*qel; % Energy gap of InAs chuang
mc = 0.023;
h=6.626070040e-34; % Planck constant, J*s
hbar=h/(2*pi); % Reduced Planck constant, J*s
m0=9.10938188E-31; % Electron mass, kg
par=Eg+hbar^2.*kxvet.^2/(2*mc*m0)*10^20;
figure,
plot(kxvet,Bands.')
hold on
plot(kxvet,par/qel,'--')

%%%
 subplot(2,1,1)
    axis on
    grid on
    hold on
    box on
    plot(mesh.xc*1e9,mesh.ecb+mesh.Eg,'b.')
    plot(mesh.xc*1e9,mesh.evb,'r.')
    xlabel('z (nm)')
    ylabel('Valence band structure, eV')
   
    subplot(2,1,2)
        axis on
        grid on
        hold on
        box on
        plot(mesh.xc*1e9,mesh.meffn,'b.')
        xlabel('z (nm)')
       ylabel('effective mass')

