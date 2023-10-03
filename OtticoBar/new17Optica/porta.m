    clear all
    close all
addpath E:\Dati\mvcsel\new10Optica
rmpath E:\Dati\mvcsel\new10Vortex    
    
    load Orta
    [Oo1,Oo2,Tort]=Teq1_modif2013Sopt(KK,lar,dos,par_grat,rr,mbvero,ifpT,segem,icriti,isca);               
           Torts{ifr}=Tort;