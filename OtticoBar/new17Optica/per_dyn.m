
  fstr=fst(:,1);
  fi=find(fstr==0);
  fstr(fi)=1;
  cstr=abs(fst(:,2));


  L_i=[];
  n_i=[];
  a_i=[];
  ics=0;
  icsc=0;
  icqw=0;
  istac=0;
  iv=1;
  clear fiC0n
  while iv<length(fstr)
%  iv
%  pausak
   for ist=1:cstr(iv)
%  ist
%  pausak
    for ist1=1:fstr(iv)
     iv1=iv-1+ist1;
     ics=ics+1;
     L_i=[L_i; dv(iv1)];
     n_i=[n_i; nv0(iv1,1)];
     a_i=[a_i; radii.a(iv1,1)];
     if (iauto(iv1,2)==-4 | iauto(iv1,2)>0)
      if istac==0
       istac=1;
      else
       icsc=icsc+1;
       fiC0n(icsc)=ics;
       istac=0;
      end
     end
%     if iauto(iv1,1)==2
     if fst(iv1,2)==-1
      icqw=icqw+1;
      istadd(icqw)=ics;
     end
     if istac==1
      icsc=icsc+1;
      fiC0n(icsc)=ics;
     end
    end
   end
   iv=iv+fstr(iv);
  end

%'controllo ficav', keyboard
  fiQ=istadd;
  fiCav=fiC0n;
  L_i=L_i/1000;


  if ifp>-1
    Ldu=[L_i*0 L_i];
    Ld=reshape(Ldu',prod(size(Ldu)),1);
    ndu=[n_i n_i];
    nd=reshape(ndu.',prod(size(Ldu)),1);
    disp(' fiQ ')
    xplo=cumsum(Ld);
    puQve=[2*fiQ-1 2*fiQ];
    puCve=[2*fiCav-1 2*fiCav];
    figure, plot(xplo,real(nd),xplo(puQve),real(nd(puQve)),'ro',...
                 xplo(puCve),real(nd(puCve)),'g.')
    keyboard
  end

if ifp~=-4
%'per_dyn', pausak
end
