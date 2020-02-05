
function reas=reconstructdi2miasd(SD,as,alpha)


nm=size(as,1);
nd=size(as,2);
sumas=sum(as');%

 for i=1:nm
 if(sumas(1,i)~=0)
    avas(i,:)=as(i,:)/sumas(1,i);
  else avas(i,:)=zeros(1,nd);  %
    end
 end
 avas=avas';
reas=zeros(nd,nm);

for i=1:nd    
     for j=1:nm
         
         for k=1:nd
               if i~=k
         reas(i,j)=reas(i,j)+alpha*SD(i,k)*avas(k,j); 
               end
           end
    end


end
      reas=as'+reas;

      