 icanc=0;
for imodi=mmvet
 raf=[nomeFs(1:end-4),num2str(iLP),num2str(imodi)];
 le=length(raf);
   ibd=0;
   ibs=0;
   ibs1=0;
   DR=dir;
   lD=length(DR);
   for kD=3:lD
    dun=getfield(DR(kD),'name');
    if length(dun)>le+4
     if isequal(dun(1:le),raf)
      eval(['hgload ',dun])
      icanc=1;
      pausak
     end
    end
   end
 if icanc==1
  vde=input(' vuoi cancellare queste figure ? [1=si] ');
  if length(vde)==1
   if vde==1
    eval(['!del ',raf,'*.fig']);
   end
  end
 end

end %imodi
