% clear 

kB=1.3806488e-23;       % Boltzmann constant (J/K)
qel=1.6021766208e-019;  % Elementary charge (C)

Px=1:mesh.nnxQW{1};       % radial nodes associated to QW
PinQW=mesh.inMQW{2}(1);   % indexes of the central QW: {2}; first column (1)

[Pnqw,Ppqw,PEfn2D,PEfp2D,PCcapn,PCcapp]=Extract_carrierQWcell(mode);

Pn3D=squeeze(MODEplot{1}.elec(:,PinQW,Px));
Pp3D=squeeze(MODEplot{1}.hole(:,PinQW,Px));
PEfn3D=squeeze(MODEplot{1}.EFn(:,PinQW,Px));
PEfp3D=squeeze(MODEplot{1}.EFp(:,PinQW,Px));

PTqw=squeeze(MODEplot{1}.Temp(:,PinQW,Px));
Pvtqw=PTqw*kB/qel;

% PCcapn=cell2mat(mode.Ccapn(2,:)');    % the first index is the QW identifier (2: central QW)
x=mesh.xgrid*1e4;
I=mode.ii_dd*1e3;
WQW=mesh.vWMQW{1};

figure
hold on,grid on,box on
plot(I,PEfn3D(:,1),'LineWidth',2)
plot(I,PEfp3D(:,1),'LineWidth',2)
plot(I,PEfn2D(:,1),'--','LineWidth',2)
plot(I,PEfp2D(:,1),'--','LineWidth',2)
title('Fermi levels')
legend('Efn 2D','Efp 2D','Efn 3D','Efp 3D')
xlabel('Current, mA'),ylabel('Energy, eV')

figure
hold on,grid on,box on
plot(I,Pn3D(:,1),'LineWidth',2)
plot(I,Pp3D(:,1),'LineWidth',2)
plot(I,Pnqw(:,1)/WQW,'LineWidth',2)
plot(I,Ppqw(:,1)/WQW,'LineWidth',2)
title('QW carriers')
legend('n 3D','p 3D','n 2D','p 2D')
xlabel('Current, mA'),ylabel('Carrier levels, cm^{-3}')

figure
hold on,grid on,box on
plot(I,PCcapn(:,1),'LineWidth',2)
plot(I,PCcapp(:,1),'LineWidth',2)
title('Capture terms')
legend('n','p')
xlabel('Current, mA'),ylabel('Capture rates, cm^{-3}\cdots^{-1}')

keyboard

figure,plot(x(Px),PEfn2D')
figure,plot(x(Px),PEfn3D')
figure,plot(x(Px),PEfn2D'-PEfn3D')
figure,plot(x(Px),(PEfn2D'-PEfn3D')./Pvtqw')
Pcap_escn=1-exp((PEfn2D'-PEfn3D')./Pvtqw');
figure,plot(x(Px),Pcap_escn),title('capture/escape')
Pfillingn=exp(-Pnqw/mode.N2{1}(1));
figure,plot(x(Px),Pfillingn),title('state filling')
Ptaun=Pn3D/mode.tausE;
figure,plot(x(Px),Ptaun),title('capture')
PCcapn_calcolato=Pfillingn.*Ptaun.*Pcap_escn';
figure,plot(x(Px),PCcapn_calcolato),title('Ccapn')
figure,plot(x(Px),PCcapn'),title('Ccapn-saved')

Pcap_escp=1-exp(-(PEfp2D'-PEfp3D')./Pvtqw');
Pfillingp=exp(-Ppqw/mode.P2{1}(1));
Ptaup=Pp3D/mode.tausH;
PCcapp_calcolato=Pfillingp.*Ptaup.*Pcap_escp';
figure,plot(x(Px),PCcapp_calcolato),title('Ccapp')
figure,plot(x(Px),PCcapp'),title('Ccapp-saved')


Prho=x;
Pdrho=[0 diff(Prho)];
Prhodrho=Prho.*Pdrho;
Prhodrho=Prhodrho(Px);

figure,plot(I,PCcapn_calcolato*Prhodrho'),title('Ccapn')
figure,plot(I,PCcapn*Prhodrho'),title('Ccapn-saved')

keyboard

[Anqw,Apqw,AEfn2D,AEfp2D]=Extract_carrierQWcell(A.mode);

An3D=squeeze(A.MODEplot{1}.elec(:,A.mesh.inMQW{2}(1),1:A.mesh.nnxQW{1}));
Ap3D=squeeze(A.MODEplot{1}.hole(:,A.mesh.inMQW{2}(1),1:A.mesh.nnxQW{1}));
AEfn3D=squeeze(A.MODEplot{1}.EFn(:,A.mesh.inMQW{2}(1),1:A.mesh.nnxQW{1}));
AEfp3D=squeeze(A.MODEplot{1}.EFp(:,A.mesh.inMQW{2}(1),1:A.mesh.nnxQW{1}));
ATqw=squeeze(A.MODEplot{1}.Temp(:,A.mesh.inMQW{2}(1),1:A.mesh.nnxQW{1}));
Avtqw=ATqw*kB/qel;

[Bnqw,Bpqw,BEfn2D,BEfp2D]=Extract_carrierQWcell(B.mode);

Bn3D=squeeze(B.MODEplot{1}.elec(:,B.mesh.inMQW{2}(1),1:B.mesh.nnxQW{1}));
Bp3D=squeeze(B.MODEplot{1}.hole(:,B.mesh.inMQW{2}(1),1:B.mesh.nnxQW{1}));
BEfn3D=squeeze(B.MODEplot{1}.EFn(:,B.mesh.inMQW{2}(1),1:B.mesh.nnxQW{1}));
BEfp3D=squeeze(B.MODEplot{1}.EFp(:,B.mesh.inMQW{2}(1),1:B.mesh.nnxQW{1}));
BTqw=squeeze(B.MODEplot{1}.Temp(:,B.mesh.inMQW{2}(1),1:B.mesh.nnxQW{1}));
Bvtqw=BTqw*kB/qel;