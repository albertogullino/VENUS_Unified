%%%%Differen.m%%%%%
global Nr_dyn
%xros=xro;
%%%%%%%%%%%%%%%%%%%%%%% non d'interesse in differen %%%%%%%%%%%%%%%
  if length(NPx)==1
   Nr_d=NPx;
  else
   Nr_d=length(NPx);
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

npfi0=NPfi;  % Zona azimutale  npfi0 = 25


  %fimaxi=2*pi;
  %fian=linspace(0,fimaxi,4*npfi0+1);

  xro=xroI; % Aggiunta per aumentare la discr. sulla coordinata radiale

  dx=diff(xro); % differenziale in r
  dxu=dx;
  dx=[dx dx(length(dx))];
  xdx=xro.*dx;

  %ifidif=1:Nr_dyn;  % da 1 a NPx1 ---> 51 o 101
  %xros=xro(ifidif); % xros == xro
  ifidif=1:length(xro);  % da 1 a NPx1 ---> 51 o 101
  xros=xro; % xros == xro

  dx1=diff(xros);
  dx1=[dx1 dx1(end)];
  xdx1=xros.*dx1;      % D'intesse per creare BD
  xdx1(end)=xdx1(end)/2;

  xvd=xro;  % compare xvd come altro nome per la discr. su r
  xu=xros; % compare xu



  ifisho=[1:npfi0+1]; % 1..26
  lfs=npfi0+1; % lfs = 26
  Nsett=4;

  if i2Ddyn==1     % 2D

     Nc=lfs; % compare Nc = 26 come valore di base
     nnu=1;

 else           % 1D
      Nc=1;
  %%%%%%%%%%%%%%%%%%%%% nn d'interesse %%%%%%%%%%%%%%%
      if i2D==3
          nnu=2;
      else
          nnu=1;
      end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end
  nnu=1;




npdif=Nc*Nr_dyn; % npdif= 26 * 51 griglia nel primo settore

  % In fian ho la discr azimutale 0 - 6.28 su 101 punti step .0628 rad

   lfi=length(fian); % 101 0 201



  if Nsett==1
   lfp=lfi-1;
   pes=ones(1,lfp);
  else                  %qui


   lfp=floor(lfi/Nsett)+1; % 26
   pes=ones(1,lfp);
   pes([1 lfp])=0.5; % riga 26 elem tutti 1 tranne 1^ e ultimo = 0.5

 end

  pesf=pes*diff(fian(1:2)); % lo moltiplico per il differenziale in fi

  if i2Ddyn==1
   lfp_N=lfp;  % 26 0 51
  else
   lfp_N=1;
  end

 if i2Ddyn==1

  xul=xu;
  pif1=1:npfi0+1;
  Npdif=length(xul);
  Nco=Npdif-2;
  Nri=length(pif1);
  pif2=[2 pif1 npfi0];
  dxf=diff(fian(1:2)); % Differenziale azimutale dxf
  nri=2:length(xul);
  Mdelf=ones(Nri,1)*xul(nri).^(-2)/dxf^2;

  dxfv=[];
  for h=1:Nri
   dxfv=[dxfv; xdx1'*pesf(h)*Nsett];
  end
  ro2D=repmat(ro',Nri,1);
%  ' qui diffe', keyboard

  imeto=1;

  if imeto==0
    Mdxr=ones(Nri,1)*(1./dx./(xul(2:Npdif)));
 %   Mdx2=ones(Nri,1)*(1./dx(1:Nco).^2);
    Mdx2=1/dx(1).^2;
  else
    xul(1)=1e-10;
    xa=xul(2)/10;
    xau=xul(Npdif)+diff(xul(Npdif-1:Npdif));

    xulp=[xul(1) xa xul(2:Npdif)];
    dxl=diff(xulp);
    Mdxr=ones(Nri,1)*(1./dxl./(xulp(2:Npdif+1)));

    dx1=[diff(xul(1:2)) diff(xul)];
 %   Mdx2=ones(Nri,1)*(1./dx1.^2);
    Mdx2=1/dx(1).^2;
    pufi=1:3;
    pufu=Npdif-3:Npdif;

  end

 else

   dxfv=xdx1'*2*pi;
   ro2D=ro';

 end %i2Ddyn

 dxfv=dxfv;

 %%%%% Mia aggiunta %%%%%%%%%%
 %fidfi=dxf.*fian(1:Nc);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%

global dxfv

 % if idiffu==0
%  imeto=1;
% end

%'differe'
%keyboard
