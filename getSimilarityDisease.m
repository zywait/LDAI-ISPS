function SD=getSimilarityDisease(interaction,SS,SSW)
%
[nd, nm] = size(interaction);
KDW=ones(nd,nd);%
for kdwi=1:nd
    for kdwj=1:nd
        if SSW(kdwi,kdwj)==0
            KDW(kdwi,kdwj)=1;
        end
        if SSW(kdwi,kdwj)==1
                KDW(kdwi,kdwj)=0;
        end
    end
end

for i=1:nd
sd(i)=norm(interaction(i,:))^2;
   end
    gamad=nd/sum(sd');
   
    %calculate Gaussian kernel for the similarity between disease: kd
    for i=1:nd
        for j=1:nd
    kd(i,j)=exp(-gamad*(norm(interaction(i,:)-interaction(j,:)))^2);
       end
    end
 
    
SD=SS.*SSW+kd.*KDW;

end
