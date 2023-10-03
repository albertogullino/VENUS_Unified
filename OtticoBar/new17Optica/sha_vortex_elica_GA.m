%'entro elica', keyboard
Par=Ps.Par;
if isfield(Par,'ScalaARC')
 ScalaARC=Par.ScalaARC;
else
 ScalaARC=0;
end

if isfield(Par,'SPP_GA')
 SPP_GA=Par.SPP_GA;
else
 SPP_GA=0;
end
% flag per scala G. Almuneau con veri settori ad indice diverso.


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
%'VER',keyboard
 if ScalaARC<=1
  nScalini=Par.Nscalini-1;
  nS2=nScalini+2;
 else
  nScalini=Par.Nscalini; 
  nS2=nScalini+1;
 end

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
  angdu=linspace(0,2*pi,nS2);
  


if iScalaContinua==0
 Ngiri=abs(Par.Ngiri);
 PhaSteps=linspace(0,2*pi,Ngiri+1);
 angduSteps=angdu;
 for k=1:Ngiri
  fi=find(angdu>=PhaSteps(k) & angdu<PhaSteps(k+1));
  angduSteps(fi)=PhaSteps(k+1);
 end
 %'qui sspi', keyboard
 if Spicchio~=0
  angSteps=angduSteps(2:end-1);
 else 
  du=fliplr(angduSteps);
  if ScalaARC<=1
   angSteps=du(3:end);
  end 
  
 end
end

%'dopo Steps', keyboard  
 if Spicchio==0
  if ScalaARC<=1
   ang_scalini=fliplr(angdu(2:end-1));
  else
    ang_scalini=fliplr(angdu);
  end
%'quo scal', keyboard
  if ScVer==1

   if ScalaARC==2   
     an_vetiVet=ang_scalini(2:end); 
     an_vetuVet=ang_scalini(1:end-1); 
   else  
    if iScalaContinua==1
     an_vetuVet=2*pi*ones(size(ang_scalini));
    else 
     an_vetuVet=angSteps;
    end
    an_vetiVet=ang_scalini;      
   end
  else
   an_vetiVet=0*ones(size(ang_scalini));
   an_vetuVet=fliplr(ang_scalini);   
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
 iScforzato=0;
 if isfield(Par,'iScforzato')==1
  iScforzato=Par.iScforzato;
  an_vetuVet=Par.angu;
  an_vetiVet=Par.angi;
  nScalini=1;
 end
 if ifp==-10
 'verifica scalino', keyboard
 end
 if ScalaARC==1
  nScalini=nScalini*2;
  lxi=nScalini;
 end  
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




if ScalaARC==1
 an_addu=an_vetiVet;
 an_addi=an_addu-an_addu(end);
 an_vetiVet=[an_vetiVet an_addi];
 an_vetuVet=[an_vetuVet an_addu];
end

if ScalaARC==2
%'fine sha_vortex_elica', keyboard
% an_addi=[an_vetiVet 0];
% an_addu=an_addi+an_vetiVet(end);
% an_vetiVet=an_addi;
% an_vetuVet=an_addu;
end

if SPP_GA==1
 nScalini=Par.Nscalini-1;
 if imag(nScalini)~=0
  if ndise==2
   nScalini=real(nScalini)-1;
  else
   nScalini=imag(nScalini)-1;
  end
 end
 PhaSteps=linspace(0,2*pi,Ngiri*(nScalini+2));
  if isfield(Par,'sfaLuna')==1
    PhaSteps=PhaSteps+Par.sfaLuna;
 %  'sfas', keyboard
  end 
 
 an_vetiVet=PhaSteps(1:end-1);
 an_vetuVet=PhaSteps(2:end);
 if Par.ScalaVerso==-1
  an_vetiVet=fliplr(an_vetiVet);
  an_vetuVet=fliplr(an_vetuVet);
 end
 lxi=nScalini;
  adis=ones(1,lxi)*Par.ri*kcav;
  pdi=zeros(1,lxi);
  ifii=zeros(1,lxi);
  bloc0=bdis(1);
  bdis=ones(1,lxi)*bloc0; 
end

%'Cont SPP', keyboard

