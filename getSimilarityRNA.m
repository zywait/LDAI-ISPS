function SM=getSimilarityRNA(interaction,FS,FSW)
%
[nd, nm] = size(interaction);
KMW=ones(nm,nm);%
for kmwi=1:nm
    for kmwj=1:nm
        if FSW(kmwi,kmwj)==0
            KMW(kmwi,kmwj)=1;
        end
        if FSW(kmwi,kmwj)==1
                KMW(kmwi,kmwj)=0;
        end
    end
end

    %calculate gamam for Gaussian kernel calculation
for i=1:nm
sm(i)=norm(interaction(:,i))^2;
   end
    gamam=nm/sum(sm');

    for i=1:nm
        for j=1:nm
    km(i,j)=exp(-gamam*(norm(interaction(:,i)-interaction(:,j)))^2);
       end
    end

SM=FS.*FSW+km.*KMW;

end
