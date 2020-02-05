
clear;
clc;

load('LncRNAdata.mat');    
 %LL:lncRNA functional similarities
 %DD:disease similarities
 alpha=0.1;%
 beta=alpha;%
 omega=0.8;%
 zeta=1-omega;
 
 nl=size(LL,1);         %nm: No. of lnRNAs
 nd=size(DD,1);         %nd: No. of diseases

ddfs=getSimilarityDisease(LD',DD,JudgeDD);
llfs=getSimilarityRNA(LD',LL,JudgeLL);

%
ldlw=reconstructmi2diasm(llfs,LD,beta);
ldlw=ldlw';
%
lddw=reconstructdi2miasd(ddfs,LD,alpha);
lddw=lddw';

ldpl=llfs*lddw; 
  for i=1:nd
  p=lddw(:,i);
     if(norm(p)~=0)
    ldpl(:,i)=ldpl(:,i)/norm(p); 
    end
  end
        
 ldpd=ddfs*ldlw;
    for i=1:nl
     p=ldlw(:,i);
      if(norm(p)~=0)
        ldpd(:,i)=ldpd(:,i)/norm(p);
      end
     end 



weight=ldpl*omega+ldpd'*zeta;

  for j=1:nd
      for i=1:nl
            
            if LD(i,j)==0  %
                 ldplnew(i,j)=ldpl(i,j);
                 ldpdnew(j,i)=ldpd(j,i);
                 prevalue(i,j)=weight(i,j);
            else   %
                 asnew=LD;  
                 asnew(i,j)=0;      
                 ddfsnew=getSimilarityDisease(asnew',DD,JudgeDD);
                 llfsnew=getSimilarityRNA(asnew',LL,JudgeLL);
                 
                 ldlwnew=reconstructmi2diasm(llfsnew,asnew,beta);
                 ldlwnew=ldlwnew';
                     % 
                 ldpdnew(j,i)=ddfsnew(j,:)*ldlwnew(:,i);
                 p=ldlwnew(:,i);
                      if(norm(p)~=0)
                     ldpdnew(j,i)=ldpdnew(j,i)/norm(p);
                      end     
                 ldpdnew(j,i)=ldpdnew(j,i);
                                                 
                    %
                 lddwnew=reconstructdi2miasd(ddfsnew,asnew,alpha);
                 lddwnew=lddwnew';
                    % 
                 ldplnew(i,j)=llfsnew(i,:)*lddwnew(:,j); 
                 
                 p=lddwnew(:,j);
                     if(norm(p)~=0)
                     ldplnew(i,j)= ldplnew(i,j)/norm(p); 
                     end
                ldplnew(i,j)=ldplnew(i,j);
                   %
                prevalue(i,j)=ldplnew(i,j)*omega+ldpdnew(j,i)*zeta;  
                
            end
        end          
    end
 

% [X_1,Y_1,tpr,aupr_1] = perfcurve(LD(:), prevalue(:),1, 'xCrit', 'reca', 'yCrit', 'prec');
[X,Y,THRE,AUC,OPTROCPT,SUBY,SUBYNAMES] = perfcurve(LD(:), prevalue(:),1);
AUC=AUC
plot(X,Y)
save result prevalue 



  
  
  
  
  
  
  
  
  
