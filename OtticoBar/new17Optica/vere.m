lambda=pi;
mm=1;
r_in=1;
rr=1;
rint=2;
r1=rint;
k0=2*pi/lambda;
kcav=k0*rr;
dos=.1;
be=1;
Idelta=.5;
KK=0;

       depe=(rint^2-rr^2)/rr^2;
       ck=-j*kcav*dos;
       Beli=ck*be;
       Mod=([Beli; -Beli]);
       
           KOt=ck*(Idelta*depe);
           KOz=0;
           P=[KOt+KOz  KOt-KOz; -(KOt-KOz)  -(KOt+KOz)]; 
           Pe=diag(Mod)+P;
           TeBW=expm(Pe);
		[V1,D1] = eig(Pe);

           Tg1=TeBW;
           
           G1=-Tg1(2,1)/Tg1(2,2);
           dd=diag(D1);
           'verifica'
           nei=abs(imag(dd(1)))/kcav/dos*rr
           ga=(r_in-r1)/(r_in+r1);
           tr=2*sqrt(r_in*r1)/(r_in+r1);
           
           Tver=1/tr*[1 -ga; -ga 1];
           Ti=inv(V1);
           Tnu=Ti/(det(Ti))^(1/4);
           fip=find(imag(dd)<0);
           fir=find(imag(dd)>0);
           fiord=[fip; fir];
           Trasn=Tnu(fiord,fiord);
           
           pu1=1:2;
           
           Vo=V1(fiord,fiord);
           Tdef=Vo*diag(exp(dd(fiord)))*inv(Vo)
           
           
            pu1=1;
            pu2=2;
            G1=-Tg1(pu2,pu1)/Tg1(pu2,pu2);    
            
            Tord=Tg1(fiord,fiord);

            T11=Tord(pu1,pu1);
            T12=Tord(pu1,pu2);
            T21=Tord(pu2,pu1);
            T22=Tord(pu2,pu2);
            
	    s12=inv(T22);
	    s21=T11-s12*T12*T21;
            s22=s12*T12;
            s11=-s12*T21;   

Tnu=Tver; 

%           ga=(r_in-r1)/(r_in+r1);
%           tr=2*sqrt(r_in*r1)/(r_in+r1);
           
%           Tver=1/tr*[1 -ga; -ga 1];
gau=-ga;
s11i=gau*exp(-j*2*dos*nei*k0);

  K=inv(eye(1)-gau*s11i)*tr;
  cdes1=[K; K*s11i]
  
  cdes2=Tver*[1; s11]
  
  %s11n=ga+tr*s11i*K
  
  %  cdes2=Tver*[1; s11n]
  

		cpiuver=zeros(1,1);
		cmenover=zeros(1,1);
		cpiuver(:,1)=[1];
		cmenover(:,1)=s11*cpiuver(:,1);
		ind=1;
	        statoin1=[cpiuver(:,ind);cmenover(:,ind)]/sqrt(r_in);;
	        
	    	statoout1=Tnu*statoin1/sqrt(r1);            
     
		ro=linspace(0,10,20);
		phi=pi/4;
		
		mpi=mm+1;
		mme=mm-1;
		Kt=rr*k0*KK;
	        fnup=cos(mpi*phi);
	    	gnup=-sin(mpi*phi);
	        fnum=cos(mme*phi);
	     	gnum=-sin(mme*phi);
		Exe=besselj(mpi,ro*Kt)*fnup+besselj(mme,ro*Kt)*fnum;
		Eye=-besselj(mpi,ro*Kt)*gnup+besselj(mme,ro*Kt)*gnum;
     		Exm=besselj(mpi,ro*Kt)*fnup-besselj(mme,ro*Kt)*fnum;
     		Eym=-besselj(mpi,ro*Kt)*gnup-besselj(mme,ro*Kt)*gnum;
     		
     		sta=statoin1;
     		SiEi=sta(1)+sta(2);
     		SiM=0;
     		
     		Exs=(SiEi*Exe+SiM*Exm);
     		Eys=(SiEi*Eye+SiM*Eym);
     		

     		sta=statoout1;
     		SiE=sta(1)+sta(2);
     		SiM=0;
     		
     		Exd=(SiE*Exe+SiM*Exm);
     		Eyd=(SiE*Eye+SiM*Eym);     		

     		figure, plot(ro,abs(Exs),ro,abs(Exd),'.'),
     		title(' Campi Ex sinistra- destra (punti) '), pausak
                
%     		figure, plot(ro,abs(Eys),ro,abs(Eyd),'.'),
%     		title(' Campi Ey sinistra- destra (punti) '), pausak
                
