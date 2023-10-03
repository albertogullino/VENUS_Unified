close all
mx=5;
my=5;
ncex=(mx-1)*2+1;
ncey=my;
Rx=.5;
Ry=.5;
dxce=2;
dyce=dxce;
Del=0;
sh_ty=0;
Tyar=1
if Tyar==-1
 Tyar=2;
end

       ncell=ncex*ncey;
       dcex=dxce;
       dcey=dyce*sqrt(3);

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

       if Tyar~=0
        ic=Tyar;
        icv=0;
        for ncc=1:ncex
         for ncr=1:ncey
          ic=ic+1;
          if is_even(ic)
           icv=icv+1;
           cceu(icv)=cced(ncr,ncc);
          end
         end
        end
%        figure,
%        plot(cce,'o'), axis equal
%        hold on, plot(cceu,'r.')
%        pausak
        cce=cceu;
       end

       fiR=find(abs(cce)<=Rad_max);
       cces=cce;
       cce=cce(fiR);

%       cce=reshape(cced,1,ncell);
       figure, plot(cce,'ro')
       grid, axis equal
       fi=linspace(0,2*pi,200);
       o=ones(size(fi))*Rad_max;
       o1=ones(size(fi))*2*dcex;
       hold on
       plot(o.*exp(j*fi),'w')
       plot(o1.*exp(j*fi),'y')
