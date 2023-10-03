colordef white 

% strName='MarkusN_TopLAqw';
% strName='btj';
% mesh.inTJ=0;

figure(105)
hold on,box on
if ~isempty(mesh.inTJ)
    % load(['standing_',strName])
%     load(['standing_',strName,'TAROCCO'])
    load standing_btjTAROCCO
    
    [AX,H1,H2]=plotyy(109+zet(end)-flip(zet), flip(nz(:,1)), 109+zet(end)-flip(zet),1.0742*flip(I0)/max(I0)*max(real(n_i)));
%     [AX,H1,H2]=plotyy(109+zet(end)-flip(zet), flip(nz(:,1)), 109+zet(end)-flip(zet),flip(I0)/241.0698*max(real(n_i)));
    
    set(AX(2),'YTick',[],'FontSize',13,'FontName','Times New Roman')
    set(AX(1),'YTick',[3 3.7],'FontSize',13,'FontName','Times New Roman')
    ylabel(AX(2),'Optical intensity, a.u.','FontSize',14,'FontName','Times New Roman')
    ylabel(AX(1),'Refractive index','FontSize',14,'FontName','Times New Roman')
    xlabel(AX(1),'z, \mum','FontSize',14,'FontName','Times New Roman')
    
    set(AX(2),'ylim',[0 4])
    set(AX(1),'ylim',[1.5 3.7])
    set(AX(2),'xlim',109+[5.6 6.6])
    set(AX(1),'xlim',109+[5.6 6.6])
    
    set(AX(1),'YColor',[0 0 1])
    set(AX(2),'YColor',[1 0 0])
    
    set(H1,'linewidth',1,'Color',[0 0 1],'linestyle','-')
    set(H2,'linewidth',1,'Color',[1 0 0],'linestyle','-')
else
%     load(['standing_',strName])
    load standing_MarkusN_TopLAqw
    [AX,H1,H2]=plotyy(109+zet(end)-flip(zet), flip(nz(:,1)), 109+zet(end)-flip(zet),flip(I0)/max(I0)*max(real(n_i)));
    
    set(AX(2),'YTick',[],'FontSize',13,'FontName','Times New Roman')
    set(AX(1),'YTick',[3 3.7],'FontSize',13,'FontName','Times New Roman')
    ylabel(AX(2),'Optical intensity, a.u.','FontSize',14,'FontName','Times New Roman')
    ylabel(AX(1),'Refractive index','FontSize',14,'FontName','Times New Roman')
    xlabel(AX(1),'z, \mum','FontSize',14,'FontName','Times New Roman')
    
    set(AX(2),'ylim',[0 4])
    set(AX(1),'ylim',[1.5 3.7])
    set(AX(2),'xlim',109+[5.6 6.6])
    set(AX(1),'xlim',109+[5.6 6.6])
    
    set(AX(1),'YColor',[0 0 1])
    set(AX(2),'YColor',[1 0 0])
    
    set(H1,'linewidth',1,'Color',[0 0 1])
    set(H2,'linewidth',1,'Color',[1 0 0])
    
    set(H1,'LineStyle','--')
    set(H2,'LineStyle','--')
    
end

grid on
