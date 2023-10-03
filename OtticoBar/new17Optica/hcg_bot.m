%'emnt', keyboard
idar=1;
Part.Refr=[1.57];
Npa1=[10];
if isfield(Par,'arc')==0
Par.arc=[20];   %spessore strato GaAs ultimo
end
%Par.primarc=[95];   %spessore prima strato GaAs 


  
alim(2)=.075; % limit of the harmonic region, comment alim(2) and dk0(2) gives co
dk0=.0015; % stepsize of k-vectors in lower k-region

alim(2)=.1; % limit of the harmonic region, comment alim(2) and dk0(2) gives co
dk0=alim(2)/Nk; % stepsize of k-vectors in lower k-region

addox=0; 



if Mis==1000
 %rad=[rad,num2str(0)];
end   
rad
%keyboard

if ipro==1
% keyboard
end 
Ps.isavetutto=(Mis+N_simf)*iSALVA  %salva tutto per indagine alvet
Ps.fiload=fiload;

if exist('i1D')
 Ps.i1D=i1D;
else
 Ps.i1D=0;
end
if Ps.i1D==1
 iBEL=100;
end
Ps.iBEL=iBEL;
 'ferma', keyboard
if exist('ifpstop')
 Ps.ifpstop=ifpstop;
else
 Ps.ifpstop=0;
end

 if ~exist('fun_wolu')
  fun_wolu='gen_woluB';
 end


if ~exist('Datt')==1
    Datt=1/2;
end

save sabard air glass rad idar Rc ifp ipro Ps Part Mis Pf Npa1 Par ifiez alim dk0 Lbuf addox ra dox fun_wolu Datt
   
%   'ver', keyboard
%bard_subud
hcg_sub_bot

load tv
global TIME   

if exist('tv')
 if length(tv)>0
  tv=tv+TIME;
 end
end
save tv tv
if ipro==2

'cont TIME', 
pausak
end

