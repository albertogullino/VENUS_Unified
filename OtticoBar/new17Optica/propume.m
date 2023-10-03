         for kr=pume
          rai=raa(kr,:);
          fao=find(rai>0);
          if length(fao)==0
           fao=0;
          end
          raa(kr,fao+1)=Ram;
          naa(kr,fao+2)=n_ext;
          shavet(kr,fao+1)=1;
         end