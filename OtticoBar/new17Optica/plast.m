clear
close all
load plast
       mx=ncex;
       my=ncey;
       ncex=(mx-1)*2+1;
       ncey=my;
       if is_even((mx-1)/2)==0
        Tyar=-1;
       else
        Tyar=1;
       end
       if Tyar==-1
        Tyar=2;
       end

       ncell=ncex*ncey;
       dcex=dxce/2;
       dcey=dyce*sqrt(3)/2;

       Rad_max=dcex*(mx-1);

       Rvetx=repmat(Rx,1,ncell);
       Rvety=repmat(Ry,1,ncell);
       Delta=repmat(Del,1,ncell);
       sh_type=repmat(sh_ty,1,ncell);

       fcex=([1:ncex]-ncex/2-.5)*dcex;
       fcey=([1:ncey]-ncey/2-.5)*dcey;
       Fcx=ones(size(fcey'))*fcex;
       Fcy=fcey'*ones(size(fcex));
       cced=Fcx+j*Fcy;

       figure, plot(cced,'ro')
       grid, axis equal
       pausak
       cceudu=ones(size(cced))*NaN;
       if Tyar~=0
        ic=Tyar;
        icv=0;
        for ncc=1:ncex
         for ncr=1:ncey
          ic=ic+1;
           if is_even(ic)
            icv=icv+1;
             cceu(icv)=cced(ncr,ncc);
             cceudu(ncr,ncc)=cced(ncr,ncc);
          end
         end
        end
        cce=cceu;
       end
       clear cceu
       figure, plot(cce,'wo')
       grid, axis equal
       hold on, plot(cceudu,'r.')
       pausak

        icdu=0;
        for ncr=1:ncey
         for ncc=1:ncex
          CEN=cceudu(ncr,ncc);
           if isnan(CEN)~=1
            icdu=icdu+1
            if icdu~=1
             icv=icv+1;
             cceu(icv)=cceudu(ncr,ncc);
             ['messo', CEN]
            else
             ['tolto', CEN]

            end
            [CEN]
            pausak
           end
           if icdu==3
            icdu=0;
           end
         end
        end
        cce=cceu;

       figure, plot(cce,'wo')
       grid, axis equal
       hold on, plot(cceudu,'r.')
       pausak

       fiR=find(abs(cce)<=Rad_max & abs(cce)~=0);
       cces=cce;
       cce=cce(fiR);

       figure, plot(cce,'r.')
       grid, axis equal
       fi=linspace(0,2*pi,200);
       o=ones(size(fi))*Rad_max;
       o1=ones(size(fi))*2*dcex;
       hold on
       plot(o.*exp(j*fi),'w')
       plot(o1.*exp(j*fi),'y')
       pausak
