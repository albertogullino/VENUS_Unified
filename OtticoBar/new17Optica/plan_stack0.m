function Tr=plan_stack0(kr,Pstack)

ie=0;
 nx=Pstack.n;
 Lx=Pstack.Li;
 rep=Pstack.rep;
 nmo=Pstack.nmo;
 kkv=Pstack.KK;
 rr=Pstack.rr;
 pol=Pstack.pol;

 be=j*kr;


% costruisco le matrici di trasmissione degli strati, per ogni kk

for ik=1:length(kkv)
    kk=kkv(ik);    
    bb=conj(sqrt(1-kk.^2));
    ZEv=(1./bb);
    ZMv=(bb);

	for in=1:length(nx)
	 Li=Lx(in);
	 ni=nx(in);
	 del=(ni^2-rr^2)/(2*rr^2)*ZEv;
	 Me=[-(bb+del) -del; del (bb+del)];
	 
	 % TE
	 Mve=Li*be*Me;
	 Tie=expm(Mve);
	 
	 % TM
	 del=(ni^2-rr^2)/(2*rr^2)*ZMv;
	 Mmt=[-(bb+del) -del; del (bb+del)];
	 delz=(1-(rr/ni)^2)/2*ZMv*(kk/bb)^2;
	 
%          P=[-(KOt+KOz)  -(KOt-KOz); 
%              (KOt-KOz)      (KOt+KOz)];	 

	 Mz=[-delz delz; -delz delz];
	 Mm=Mmt+Mz;
	 Mvm=Li*be*Mm;
	 Tim=expm(Mvm); 
	
	 Tme(:,:,in)=Tie;
	 Tmm(:,:,in)=Tim;
%' interno strati', keyboard	 
	end

%' interno k stack0', keyboard

  irep=1;
  Te=eye(2);
  Tm=eye(2);
  while irep<=length(nx)
%   'inizio while',    irep
   if rep(irep)==0
    Te=Tme(:,:,irep)*Te;
    Tm=Tmm(:,:,irep)*Tm;
    irep=irep+1;
   else
     TeM=eye(2);
     TmM=eye(2);
     NPairs=rep(irep);
     while rep(irep)==NPairs
      TeM=Tme(:,:,irep)*TeM;
      TmM=Tmm(:,:,irep)*TmM;     
      irep=irep+1;
%      irep, pausak
      if irep>length(rep)
       break
      end
     end
    Te=TeM^NPairs*Te;
    Tm=TmM^NPairs*Tm;
    TeMi=TeM^NPairs;
    TmMi=TmM^NPairs;    
    TeMi=TeM;
    TmMi=TmM;    
%    'qui mir', keyboard
   end
   
  end
%' fine k', keyboard
 Vepp(ik,1)=Te(1,1);
 Verr(ik,1)=Te(2,2);
 Vepr(ik,1)=Te(1,2);
 Verp(ik,1)=Te(2,1);
 
 Vmpp(ik,1)=Tm(1,1);
 Vmrr(ik,1)=Tm(2,2);
 Vmpr(ik,1)=Tm(1,2);
 Vmrp(ik,1)=Tm(2,1);
 if exist('TeMi')
 ie=1;
 Veppm(ik,1)=TeMi(1,1);
 Verrm(ik,1)=TeMi(2,2);
 Veprm(ik,1)=TeMi(1,2);
 Verpm(ik,1)=TeMi(2,1);
 
 Vmppm(ik,1)=TmMi(1,1);
 Vmrrm(ik,1)=TmMi(2,2);
 Vmprm(ik,1)=TmMi(1,2);
 Vmrpm(ik,1)=TmMi(2,1); 
 end
 
end  %ik
VeppM=repmat(Vepp,nmo,1);
VerrM=repmat(Verr,nmo,1);
VeprM=repmat(Vepr,nmo,1);
VerpM=repmat(Verp,nmo,1);

VmppM=repmat(Vmpp,nmo,1);
VmrrM=repmat(Vmrr,nmo,1);
VmprM=repmat(Vmpr,nmo,1);
VmrpM=repmat(Vmrp,nmo,1);

 Mpp=diag([VeppM; VmppM]);
 Mrr=diag([VerrM; VmrrM]);
 Mrp=diag([VerpM; VmrpM]);
 Mpr=diag([VeprM; VmprM]);
 Z=zeros(size(Mpp));
 Pus0=Pstack.Pus0;

if pol==0

 M11 =[Mpp Z; ... 
 	 Z Mpp]; 
 M12 =[Mpr Z; ... 
 	  Z Mpr];  	 
 M22 =[Mrr Z; ... 
 	 Z Mrr]; 
 M21 =[Mrp Z; ... 
 	  Z Mrp];  	  	  
% Tr=[Mpp Z Mpr Z; ...
%       Z Mpp Z Mpr; ...
%       Mrp Z Mrr  Z; ...
%       Z Mrp Z Mrr]; 
 Tr=[ M11(Pus0,Pus0) M12(Pus0,Pus0); ...
        M21(Pus0,Pus0) M22(Pus0,Pus0)];
       
else
% Vd=[VeppM; VmppM; VerrM; VmrrM];
% Vdp=[VeprM; VmprM];
% Vdm=[VerpM; VmrpM];
%
% le=length(Vdp);
% Tr=diag(Vd)+diag(Vdp,le)+diag(Vdm,-le);
% 
 Tr=[Mpp(Pus0,Pus0) Mpr(Pus0,Pus0); Mrp(Pus0,Pus0) Mrr(Pus0,Pus0)];
 
end

%' ferma stack0', keyboard


% verifiche per estrarre i vari termini: non utilizzato ora
if ie==0
return
end

VeppMi=repmat(Veppm,nmo,1);
VerrMi=repmat(Verrm,nmo,1);
VeprMi=repmat(Veprm,nmo,1);
VerpMi=repmat(Verpm,nmo,1);

VmppMi=repmat(Vmppm,nmo,1);
VmrrMi=repmat(Vmrrm,nmo,1);
VmprMi=repmat(Vmprm,nmo,1);
VmrpMi=repmat(Vmrpm,nmo,1);

Vdm=[VeppMi; VmppMi; VerrMi; VmrrMi];
Vdpm=[VeprMi; VmprMi];
Vdmm=[VerpMi; VmrpMi];

le=length(Vdpm);
TrMi=diag(Vdm)+diag(Vdpm,le)+diag(Vdmm,-le);

