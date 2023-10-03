if mode.EqPermutazione==0
    %
    nn=mesh.nn;                   % -> elec eqs.
    pp=nn+mode.nflg*nn;           % -> hole eqs.
    tt=pp+mode.pflg*nn;           % -> trap eqs.
    qq=tt+mode.tflg*nl*nn;        % -> circuit eqs.
    rr=qq+nm;                     % -> current eqs.
    vv=rr+nm;                     % -> n2D eqs.
    ww=vv+length(mesh.xgrid)*NQW;               % -> p2D eqs.
    ss=ww+length(mesh.xgrid)*NQW;               % -> optical rate eqs.
else
    nn=mesh.nn;                             % -> elec eqs.
    pp=nn+mode.nflg*nn;                     % -> hole eqs.
    tt=pp+mode.pflg*nn;                     % -> trap eqs.
    vv=tt+mode.tflg*nl*nn;                  % -> n2D eqs
    ww=vv+length(mesh.xgrid)*NQW;           % -> p2D eqs
    ss=ww+length(mesh.xgrid)*NQW;           % -> optical rate eqs.
    qq=ss+mode.nmodes;
    rr=qq+nm;
end
neq=nn+mode.nflg*nn+mode.pflg*nn+mode.tflg*nl*nn+2*nm+mode.oflg*not(mode.firstrun)*2*NQW*length(mesh.xgrid)+mode.oflg*mode.nmodes;
tvetAC = sparse(neq,1);
tvetAC(qq+1) = 1;

JmatDC = Jmat0;

Y_ss = zeros(1,length(mode.fvet));
G_ss = zeros(1,length(mode.fvet));
C_ss = zeros(1,length(mode.fvet));
Pst_ss = zeros(mode.nmodes,length(mode.fvet));


% for iomega = 1:length(mode.fvet)
parfor iomega = 1:length(mode.fvet)
%     iomega = mode.fvet(1);

    omega = 2*pi*mode.fvet(iomega);
    JmatAC = JmatDC + 1i.*omega.*Jmat1;
    
    delta = JmatAC\tvetAC;
    
    Y=delta(rr+1); % small-signal admittance, ohm/cm^2
    
    Y_ss(iomega)=Y;
    G_ss(iomega)=real(Y); % small-signal conductance, ohm/cm^2
    C_ss(iomega)=imag(Y/omega);  % small-signal capacitance, ohm/cm^2
%     Z_ss(iomega)=1./Y;
%     R_ss(iomega)=real(Z_ss);
    
    % Parasitic Parameters
    %     Rm = 80;
    %     Cp = 900e-15;
    %     Radd = 0;
    %
    % NUSOD 2020
    Rm = 4;
    Cp = 15.5e-12;
    fRC = 1/(2*pi*Rm*Cp)*1e-9
    
    Ypad = 1/((1/Y)+Rm);
    Ypar = Ypad + 1i*omega*Cp;
    
    Pst_ss(:,iomega) = delta(ss+[1:mode.nmodes]);
    
    disp(['Small-signal analysis ',num2str(iomega),' of ',num2str(length(mode.fvet)),' executed'])
    
end
%
% AM
% AM = abs(sum(Pst_ss)./Y_ss);
AM = abs(sum(Pst_ss)./Y_ss).^2;
% AM = AM./AM(1);
% S11
% Z0 = 50;
% S11 = (Z*mode.Area-Z0)./(Z*mode.Area+Z0);
%

% return
if IPLOT==1
    figure
    set(gcf,'Position',[379 522 1079 420])
    subplot(1,3,1)
    semilogx(mode.fvet/1e9,G_ss,'.-')
    grid on
    hold on
    xlabel('Frequency, GHz')
    ylabel('Differential conductance, S')
    subplot(1,3,2)
    semilogx(mode.fvet/1e9,C_ss,'.-')
    grid on
    hold on
    xlabel('Frequency, GHz')
    ylabel('Differential capacitance, F')
    subplot(1,3,3)
    semilogx(mode.fvet/1e9,10*log10(AM/AM(1)),'.-')
    grid on
    hold on
    xlabel('Frequency, GHz')
    ylabel('AM response, normalized')
    title(['Current: ',num2str(Cf),' mA'])
end
% keyboard