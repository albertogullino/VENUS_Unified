

       ncell=ncex*ncey;
       dcex=2*Rx+dxce;
       dcey=2*Ry+dyce;
       Rvetx=repmat(Rx,1,ncell);
       Rvety=repmat(Ry,1,ncell);
       Delta=repmat(Del,1,ncell);
       sh_type=repmat(sh_ty,1,ncell);

       fcex=([1:ncex]-ncex/2-.5)*dcex;
       fcey=([1:ncey]-ncey/2-.5)*dcey;
       Fcx=ones(size(fcey'))*fcex;
       Fcy=fcey'*ones(size(fcex));
       cced=Fcx+j*Fcy;

       cce=reshape(cced,1,ncell);

       ic=Tyar;
       icv=0;
       for ncc=1:ncex
        for ncr=1:ncey
         ic=ic+1;
         if is_even(ic)
          icv=icv+1;
          cceu(icv)=cced(ncc,ncr);
         end
        end
       end
       figure,
       plot(cce,'o'), axis equal
       hold on, plot(cceu,'r.')
