%
%  Versione 1.1 10-Febbraio-2006 di Oscar A. Peverini
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Sintassi:  [xvet,wvet]=quadad(tipo,a,N,ALPHA,BETA)
%
% Ingressi:
%         tipo-----stringa che identifica il metodo di quadratura:
%                  'lague'---Gauss-Laguerre tra 0 ed Inf;
%                  'legen'---Gauss-Legendre tra -a e a;
%                  'legsh'---Gauss-Legendre shifted, ossia tra 0 e a;
%                  'cheby'---Gauss-Chebyshev di 1°specie tra -1 e 1;
%                  'jacob'---Gauss-Jacobi con peso (1-x)^ALPHA*(1+x)^BETA con parametri ALPHA e BETA>-1 tra -1 e 1;
%         a--------estremo dell' intervallo di integrazione [-a,a] di integrazione;
%                  nel caso di tipo='lague' e 'cheby' occorre porre a=1;
%         N--------numero di nodi della formula di quadratura;
%       
% Uscite:
%         xvet-----vettore contenente gli zeri della formula di quadratura;
%         wvet-----vettore contenente i pesi della formula di quadratura;
% Esempio:
%         » N=10;a=2.3;[xvet,wvet]=quadad('legen',a,N);ris1=sum(wvet.*cos(xvet))
%         ris1=1.4914
%
%         » ris2=quad8('cos',-a,a)
%         ris2=1.4914
%         »
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         

function[xvet,wvet]=quadad(tipo,a,N,ALPHA,BETA);


if(tipo=='jacob' & ALPHA==0 & BETA==0)
    tipo='legen';
end

m=0:1:N-2;
if tipo=='lague'
   alpha=-(m+1);
elseif tipo=='legen'
   alpha=(m+1)./(2*m+1).*sqrt((2*m+1)./(2*m+3));   
elseif tipo=='legsh'
   alpha=(m+1)./2./(2*m+1).*sqrt((2*m+1)./(2*m+3));   
elseif tipo=='cheby';
  m=1;     
elseif tipo=='jacob'
    hm=2^(ALPHA+BETA+1)./(2*m+ALPHA+BETA+1).*gamma(m+ALPHA+1).*gamma(m+BETA+1)./gamma(m+1)./gamma(m+ALPHA+BETA+1);
    hmpiu1=2^(ALPHA+BETA+1)./(2*(m+1)+ALPHA+BETA+1).*gamma((m+1)+ALPHA+1).*gamma((m+1)+BETA+1)./gamma((m+1)+1)./gamma((m+1)+ALPHA+BETA+1);
    a1m=2*(m+1).*(m+ALPHA+BETA+1).*(2*m+ALPHA+BETA);
    a2m=(2*m+ALPHA+BETA+1)*(ALPHA^2-BETA^2);
    a3m=gamma((2*m+ALPHA+BETA)+3)./gamma((2*m+ALPHA+BETA));
    alpha=a1m.*sqrt(hmpiu1)./(a3m.*sqrt(hm));
%    alpha=a1m./a3m;
end
m=0:1:N-1;
if tipo=='lague'
   beta=(2*m+1);
elseif tipo=='legen'
   beta=0*m;
elseif tipo=='legsh'
   beta=1/2*ones(1,length(m));  
elseif tipo=='cheby'
   m=1;
elseif tipo=='jacob'
    hm=2^(ALPHA+BETA+1)./(2*m+ALPHA+BETA+1).*gamma(m+ALPHA+1).*gamma(m+BETA+1)./gamma(m+1)./gamma(m+ALPHA+BETA+1);
    hmpiu1=2^(ALPHA+BETA+1)./(2*(m+1)+ALPHA+BETA+1).*gamma((m+1)+ALPHA+1).*gamma((m+1)+BETA+1)./gamma((m+1)+1)./gamma((m+1)+ALPHA+BETA+1);
    a1m=2*(m+1).*(m+ALPHA+BETA+1).*(2*m+ALPHA+BETA);
    a2m=(2*m+ALPHA+BETA+1)*(ALPHA^2-BETA^2);
    a3m=gamma((2*m+ALPHA+BETA)+3)./gamma((2*m+ALPHA+BETA));
    beta=-a2m./a3m;
end;

if(tipo=='lague'|tipo=='legen'|tipo=='legsh'|tipo=='jacob')
    T=diag(alpha,1)+diag(beta,0)+diag(alpha,-1);
    [U,D]=eig(T);
    xvet=diag(D).';
    wvet=U(1,:).^2;
    if(tipo=='jacob')
        wvet=wvet*hm(1)/sum(wvet);
    end
end
if tipo=='legen'
    wvet=wvet*2;
elseif tipo=='cheby'
    wvet=pi/N*ones(1,N);
    n=1:N;
    xvet=cos((2*n-1)/2/N*pi);
elseif tipo=='chemo'
    [xvet0,wvet0]=quadad('legen',1,N,0,0);
    [xvet0,I]=sort(xvet0);
    wvet0=wvet0(I);
    k=1:N;
    thk=(2*k-1)*pi/2/N;
    xk=cos(thk);
    wk=zeros(1,N);
    w=(xvet0+1).*(xvet0-1);
    for k=1:N
        Lk=ones(1,length(xvet0));
        for i=1:N
            if(i~=k)
                Lk=Lk.*(xvet0-xk(i))./(xk(k)-xk(i));
            end
        end

        wk(k)=(Lk)*wvet0.';

        %         figure(1)
        %         hPl=plot(xk,zeros(1,N),'ro',xvet0,Lk,'b:',xvet0,w,'g-');
        %         set(hPl,'LineWidth',1.2);
        %         grid on
        %         title(['k = ',num2str(k)]);
        %         drawnow
        %         pause;%(0.01)
    end
    xvet=xk;
    wvet=wk;
end

xvet=xvet*a;
wvet=wvet*a;