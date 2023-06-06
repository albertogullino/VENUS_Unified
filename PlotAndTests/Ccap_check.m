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

figure,plot(mesh.xgrid(Px)*1e4,PEfn2D')
figure,plot(mesh.xgrid(Px)*1e4,PEfn3D')
figure,plot(mesh.xgrid(Px)*1e4,PEfn2D'-PEfn3D')
figure,plot(mesh.xgrid(Px)*1e4,(PEfn2D'-PEfn3D')./Pvtqw')
Pcap_escn=1-exp((PEfn2D'-PEfn3D')./Pvtqw');
figure,plot(mesh.xgrid(Px)*1e4,Pcap_escn),title('capture/escape')
Pfillingn=exp(-Pnqw/mode.N2{1}(1));
figure,plot(mesh.xgrid(Px)*1e4,Pfillingn),title('state filling')
Ptaun=Pn3D/mode.tausE;
figure,plot(mesh.xgrid(Px)*1e4,Ptaun),title('capture')
PCcapn_calcolato=Pfillingn.*Ptaun.*Pcap_escn';
figure,plot(mesh.xgrid(Px)*1e4,PCcapn_calcolato),title('Ccapn')
figure,plot(mesh.xgrid(Px)*1e4,PCcapn'),title('Ccapn-saved')

Pcap_escp=1-exp(-(PEfp2D'-PEfp3D')./Pvtqw');
Pfillingp=exp(-Ppqw/mode.P2{1}(1));
Ptaup=Pp3D/mode.tausH;
PCcapp_calcolato=Pfillingp.*Ptaup.*Pcap_escp';
figure,plot(mesh.xgrid(Px)*1e4,PCcapp_calcolato),title('Ccapp')
figure,plot(mesh.xgrid(Px)*1e4,PCcapp'),title('Ccapp-saved')


Prho=mesh.xgrid*1e4;
Pdrho=[0 diff(rho)];
Prhodrho=Prho.*Pdrho;
Prhodrho=Prhodrho(Px);

figure,plot(mode.ii_dd*1e3,PCcapn_calcolato*Prhodrho'),title('Ccapn')
figure,plot(mode.ii_dd*1e3,PCcapn*Prhodrho'),title('Ccapn-saved')


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