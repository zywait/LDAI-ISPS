
function reas=reconstructmi2diasm(SM,as,alpha)


nm=size(as,1);
nd=size(as,2);

sumas=sum(as);

 for i=1:nd
 if(sumas(1,i)~=0)
    avas(:,i)=as(:,i)/sumas(1,i);
  else avas(:,i)=zeros(nm,1); 
    end
 end



 
reas=zeros(nm,nd);

for i=1:nm    
     for j=1:nd
         
         for k=1:nm
               if i~=k
          reas(i,j)=reas(i,j)+alpha*SM(i,k)*avas(k,j); 
               end
           end
    end


end

    reas=as+reas;

      