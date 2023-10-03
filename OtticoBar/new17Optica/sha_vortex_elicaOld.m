
Par=Ps.Par;

%'entro elica', keyboard

if isfield(Par,'SPP')==1
 SPP=Par.SPP;
else 
 SPP=0;
end

if isfield(Par,'iScalaContinua')==1
  iScalaContinua=Par.iScalaContinua;
else
  iScalaContinua=1;
end

isimm=0;

if SPP==0
 sha_vortex
else
 nScalini=Par.Nscalini-1;
 if isfield(Par,'ScalaVerso')==0
  ScVer=1;
 else
  ScVer=Par.ScalaVerso;
 end
 if Par.Spicchio~=0 
  Spicchio=Par.Spicchio;
  lxi=1;
 else 
  lxi=nScalini;
  Spicchio=0;
 end
  angdu=linspace(0,2*pi,nScalini+2);
  


if iScalaContinua==0
 Ngiri=abs(Par.Ngiri);
 PhaSteps=linspace(0,2*pi,Ngiri+1);
 angduSteps=angdu;
 for k=1:Ngiri
  fi=find(angdu>=PhaSteps(k) & angdu<PhaSteps(k+1));
  angduSteps(fi)=PhaSteps(k+1);
 end
 if Spicchio==1
  angSteps=angduSteps(2:end-1);
 else 
  du=fliplr(angduSteps);
  angSteps=du(3:end);
 end
end

%'dopo Steps', keyboard  
 if Spicchio==0
  ang_scalini=fliplr(angdu(2:end-1));
%'quo scal', keyboard
  if ScVer==1
   if iScalaContinua==1
    an_vetuVet=2*pi*ones(size(ang_scalini));
   else 
    an_vetuVet=angSteps;
   end
   an_vetiVet=ang_scalini; 
  else
   an_vetiVet=0*ones(size(ang_scalini));
   an_vetuVet=ang_scalini;   
  end
  diAn=an_vetiVet-an_vetuVet;
  fi_0=find(diAn~=0);
  an_vetiVet=an_vetiVet(fi_0);
  an_vetuVet=an_vetuVet(fi_0);
  nScalini=length(fi_0);
  lxi=nScalini;
 else
  ang_scalini=angdu(1:end);
  an_vetuVet=2*pi*ones(size(ang_scalini));
  an_vetiVet=ang_scalini(Spicchio);
  an_vetuVet=ang_scalini(Spicchio+1);
  nScalini=1;
 end
 
% 'verifica scalino', keyboard
  adis=zeros(1,lxi);
  pdi=zeros(1,lxi);
  ifii=zeros(1,lxi);
  bloc0=bdis;
  bdis=ones(1,lxi)*bloc0;
  ifalso=-1
%  keyboard
end

if isimm==1
 dfi=an_vetiVet(end)*4;
 an_vetiVet= an_vetiVet+dfi;
 an_vetuVet= an_vetuVet+dfi;
end