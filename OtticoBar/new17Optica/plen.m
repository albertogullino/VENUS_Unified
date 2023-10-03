%load le
close all
clear
%load levs
load le

Rm_rel=0
   [dilu,ailu,nilu,filu]=lens_suba(duth,ha,ra,Ndisc,dunr,UD,-11,Rel,Nrel,sNpa*Npair,Rflat,Rm_rel,Np_ad);
%   [dilu,ailu,nilu,filu]=lens_subc(duth,ha,ra,Ndisc,dunr,UD,-11,Rel,Nrel,sNpa*Npair,Rflat,Rm_rel);

