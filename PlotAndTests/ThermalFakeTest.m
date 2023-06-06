clear 
close all
clc

Ivet=linspace(0,13,100) ; 
load('Tfit')

for indI=1:100
     Current=Ivet(indI) ; 
     indM=find(Current>=I) ;
     indM=indM(end) ; 
     I1=I(indM) ; 
     I2=I(indM+1) ; 
     DeltaI=I2-I1 ;
     DeltaT=squeeze( DT(indM+1,:,:)-DT(indM,:,:) ) ; 
     slope=DeltaT./DeltaI ;
     DeltaIfit=Current-I1 ;
     DeltaTold=squeeze(DT(indM,:,:))+slope.*DeltaIfit ;
     Tprec=293+DeltaTold ;  
     DeltaTmax(indI)=max(max(DeltaTold)) ; 
     
    
    
    
end