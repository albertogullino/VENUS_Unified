%  'prima deter_gamma in Ass_lente: vedi Gac  ', keyboard
ireturn=0;
if exist('gra_le')
 if isfield(gra_le,'pos')==0
   gra_le.Gac=0;
   gra_le.ru=0;
   Gac=0;
   ireturn=1;
 else 
  if gra_le.pos==1
   clear global igraef_new
   igraef_new=0;
%   ' modificato _new', keyboard
  end
 end 
end  
%' prima di det_gamma', keyboard
 deter_gamma
 if ifp==-10
  'dopo deter_gamma in ass_parn: vedi Gac  '
  pausak
 end

  if Npair>0
%   [dilu,ailu,nilu,ipmir]=lens_sub(duth,ha,ra,Ndisc,dunr,UD,ifp,Rel,Nrel,Npair);
   [dilu,ailu,nilu,filu]=lens_sub(duth,ha,ra,Ndisc,dunr,UD,ifp,Rel,Nrel,Npair,Rflat,Rm_rel,gra_le);
   no_pair=1;
   no_rep=0;
%   pu_fst=[ones(size(dilu))*no_rep ones(size(dilu))*no_pair];
   pu_fst=filu;
  else
   per0=sum(duth(1:2))/1000;
   Npe=floor(2*ha/per0)+1;
   if Npair-floor(Npair)>0
    Npadu=Npe+.5;
   else
    Npadu=Npe;
   end
%   ' prima lands', keyboard
%   [dilu,ailu,nilu,ipmir]=lens_sub(duth,ha,ra,Ndisc,dunr,UD,ifp,Rel,Nrel,Npadu);
   Npadu=Npair;
   if isfield(PAle{klen},'Np_ad')==1
%        'problema 2 Np_ad in ass_parn', keyboard   
     Np_ad=PAle{klen}.Np_ad;
     if length(Np_ad)>0
       if Np_ad<0
         Npadu(2)=Np_ad;        
       else
        'problema 1 Np_ad in ass_parn', keyboard
       end
      else
%        'problema 2 Np_ad in ass_parn', keyboard
      end
   end

%        'NP_ad in ass_parn', keyboard
%  'prima kessub 2', keyboard

   [dilu,ailu,nilu,filu,ngvero,gra_le]=lens_sub(duth,ha,ra,Ndisc,dunr,UD,ifp,Rel,Nrel,Npadu,Rflat,Rm_rel,gra_le);
%   no_pair=floor(abs(Npair))-Npe+1;
%   no_rep=length(ipmir);
%   pu_fst=[ones(size(dilu))*0 ones(size(dilu))*1];
%   pu_fst(ipmir,1)=no_rep;
%   pu_fst(ipmir,2)=no_pair;
    pu_fst=filu;
  end
  %'dopo kessub', keyboard
  
  Dilu{klen}=dilu;
  Nilu{klen}=nilu;
  Ailu{klen}=ailu;
  Filu{klen}=filu;
  Raxv(klen)=Rax;
  ShMv(klen)=ShaM;
