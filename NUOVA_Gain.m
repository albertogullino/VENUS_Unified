global N2lin P2lin N2i P2i T2 L2 indGan
global G0 dgE0 dgH0 rsp0 drspE0 drspH0 Rsp0 dRspE0 dRspH0 DeltaN0

tic
meshT=mesh;
meshT.T=mode.T0+(mesh.Tgain-mode.T0)*mode.C_TempGain;
[N2lin,P2lin,N2i,P2i,T2,L2,indGan]=vet_mat(mode,meshT);

%'vedo NUOVA', keyboard
% 15/12/2017
% Prima di oggi, facevamo tutto vettoriale (avevamo che L2 era un
% vettore/matrice Nmodi x Nelementi) perché non avevamo capito benissimo
% come usare interpn. Usando una/due lunghezze d'onda, conviene fare come
% ora. Vedi NUOVA_Gain_OLD_BACKUP per come si faceva prima.
% In questo modo è molto più veloce.
% Il codice quindi è un po' raffazzonato perché è stato modificato meno
% possibile, ma si può mettere più elegante (evitare di mettere L2 come
% vettore lungo / matrice). Allo stesso modo, non ha più senso che N2lin,
% P2lin, T2 siano matrici: è sufficiente che siano sempre vettori di
% lunghezza Nelementi.

G0=zeros(size(N2lin));
dgE0=G0;
dgH0=G0;
rsp0=G0;
drspE0=G0;
drspH0=G0;
DeltaN0=G0;
for indlam=1:size(L2,1)
    np_E=N2lin-N2i;
    np_H=P2lin-P2i;
    %'ver qui NUOVA', keyboard
    [Rsp0,dRspE0,dRspH0]=f_InterpRsp4D(np_E,np_H,T2);
%     tic
    [vG0,vdgE0,vdgH0,vrsp0,vdrspE0,vdrspH0,vDeltaN0] =...
        f_InterpGain4D(N2lin(1,:),P2lin(1,:),mode.Deltalam+L2(indlam,1),T2(1,:));
%     toc
    for ind_diag=1:length(N2lin)
        G0(indlam,ind_diag)=vG0(ind_diag,ind_diag,1,ind_diag);
        dgE0(indlam,ind_diag)=vdgE0(ind_diag,ind_diag,1,ind_diag);
        dgH0(indlam,ind_diag)=vdgH0(ind_diag,ind_diag,1,ind_diag);
        rsp0(indlam,ind_diag)=vrsp0(ind_diag,ind_diag,1,ind_diag);
        drspE0(indlam,ind_diag)=vdrspE0(ind_diag,ind_diag,1,ind_diag);
        drspH0(indlam,ind_diag)=vdrspH0(ind_diag,ind_diag,1,ind_diag);
        DeltaN0(indlam,ind_diag)=vDeltaN0(ind_diag,ind_diag,1,ind_diag);
    end

end



if length(find(isnan(G0)==1))>0
 'G0 anan', keyboard
end

if v0_dd>1
    NQW=mesh.NMQW;
    sDe=size(DeltaN0);
    DeDe=mean(DeltaN0,1);   % no effect in 1D
    Deld=reshape(DeDe,sDe(2)/NQW,NQW);
    mode.DeltaN=mean(Deld,2);
%     if mode.quasi1D==1
%         mode.DeltaN=DeltaN0;
%     end
else
    mode.DeltaN=0;
end

% toc
 %' in nuova Gain', keyboard


