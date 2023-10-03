% neff: indice efficace TE neff(1) e TM neff(2)
% err: errore sul secondo termine della matrice catena
% G: coeff. riflessione TE e TM corrispondenti

function [neff,err,G,nefb]=neffic(tetai,r_in,r_out,n1,n2,d1i,d2i,dx,lambda0,NModi);

    [Te,Tm,Ge,Gm,du,du,du,du,nefb]=Orta_tracar(tetai,r_in,r_out,n1,n2,d1i,d2i,dx,lambda0,NModi,3);
% [Te,Tm,Ge,Gm]=Orta_tracar(tetai,r_in,r_out,n1,n2,d1i,d2i,dx,lambda0,NModi,itetmt,Strut);

%  [Te,Tm,Ge,Gm]=Orta_tras(tetai,r_in,r_out,n1,n2,d1i,d2i,dx,lambda0,NModi,itetmt);
%  [Te,Tm,Ge,Gm,du,du,du,du,nef]=bast_tras(tetai,r_in,r_out,n1,n2,d1i,d2i,dx,lambda0,NModi,isto,itetmt);
    G(1)=Ge;
    G(2)=Gm;
    rr=1;
    Mu=[1 1; r_out/rr -r_out/rr]/sqrt(r_out); 
    Mi=0.5*sqrt(r_in)*[1 rr/r_in; 1 -rr/r_in];    
    KKie=(Mu*Te*Mi);
    KKim=(Mu*Tm*Mi);
    k0=2*pi/lambda0;
    neff(1)=acos(KKie(1,1))/(dx*k0);  
    if imag(neff(1))>0
     neff(1)=real(neff(1))-j*imag(neff(1));
    end
    if abs(imag(neff(1)))>1e6*abs(real(neff(1)))
      neff(1)=-j*abs(imag(neff(1)));
    end
    neff(2)=acos(KKim(1,1))/(dx*k0);  
    if imag(neff(2))>0
     neff(2)=real(neff(2))-j*imag(neff(2));
    end    
    if abs(imag(neff(2)))>1e6*abs(real(neff(2)))
      neff(2)=-j*abs(imag(neff(2)));
    end
    fas=dx*k0*neff;
    mi=rr./neff;
    K12=-j*sin(fas).*mi;
    K12
    Kv12=[KKie(1,2) KKim(1,2)]
    err=ones(size(K12))-abs(K12./Kv12); 
    if length(find(err>.5))>0
     'indice efficace non valido'
     fin=find(err>.5);
     neff(fin)=0;
    end
%'in neffic',     keyboard