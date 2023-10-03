    if length(xroJ)==1
     Nx=xroJ;
     plomax=2.5*max(max(axtot));
     plomay=2.5*max(max(aytot));
     ploma=max([plomax plomay]);
     xvero=linspace(0,ploma,Nx);
     x=xvero;
    else
     xvero=xroJ;
     x=xvero;
    end

%' xvero  ', keyboard

%    disp(' camful'), keyboard

    if length(lfi_inp)==0 | length(fimaxi)==0
     lfi_inp=51;
     fimaxi=2*pi;
    end

    Nx=length(x);
    npk=length(KK);


  lbv=length(mbv);
  fian=linspace(0,fimaxi,lfi_inp(1)+1);
%  fian=linspace(0,fimaxi,lfi_inp(1));
  if fimaxi<2*pi
   fian0=[fian(1:length(fian)-1) fian+pi];
  else
   fian0=fian;
  end

  fm0=cos(fian0);
  gm0=sin(fian0);

    xdx=[0 diff(xvero)].*xvero;
    defid=diff(fian0);
    defi=[defid defid(end)]';
    defi([1 end])=defi([1 end])/2;
    XP=xvero'*fm0;
    YP=xvero'*gm0;


    iLP=iLPr;
%    iLP=iLP1;
pola=ipolar;
%' modangle', keyboard
ipost=0
    modangle
    mod_ff
