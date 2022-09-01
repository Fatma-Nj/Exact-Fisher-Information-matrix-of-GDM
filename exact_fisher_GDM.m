function fisher= exact_fisher_GDM(X,K,pi,theta,p,tau)
%% calculating the exact fisher information matrix of the generalized dirichlet multinomial
[N, D] = size(X);
fisher=cell(K,D-1);

posterior=posterior_GDM(X,K,pi,theta,p,tau);
for j=1: K
  for d=1:D-1        
        fisher{j,d}=zeros(2,2);fish1=0;fish2=0;fish3=0;
     for i=1:N    
        Y=sum(X(i,d:D));
        pdf_X=pdf_gdm(X(i,d),K,pi(j,d),theta(j,d));
        Y_l=sum(X(i,d+1:D));
        pdf_Y=pdf_gdm(Y,K,pi(j,d),theta(j,d));
        pdf_Yl=pdf_gdm(Y_l,K,pi(j,d),theta(j,d));
        m=sum(X(i,:),2);sum_X1=0; sum_X2=0; sum_X3=0; sum_Y1=0;sum_Y2=0;sum_Y3=0;sum_Y4=0;
       for r=1:m
        if (X(i,d) >= r)
         sum_X1 = sum_X1 + pdf_X(:,j)./(pi(j,d) + (r-1).* theta(j,d)).^2;
         sum_X2 = sum_X2 + ((r-1).* pdf_X(:,j))./(pi(j,d) + (r-1).* theta(j,d)).^2;
         sum_X3 = sum_X3 + ((r-1).^2.* pdf_X(:,j))./(pi(j,d) + (r-1).* theta(j,d)).^2;
        end

        if (Y >=r)
         sum_Y1 = sum_Y1+ ((r-1).^2.* pdf_Y(:,j))./(1 + (r-1).* theta(j,d)).^2;
        end
      
        if (Y_l >=r)
         sum_Y2 = sum_Y2+ pdf_Yl(:,j)./(1-pi(j,d) + (r-1).* theta(j,d)).^2;    
         sum_Y3 = sum_Y3+ ((r-1).* pdf_Yl(:,j))./(1-pi(j,d) + (r-1).* theta(j,d)).^2;   
         sum_Y4 = sum_Y4+ ( (r-1).^2 .* pdf_Yl(:,j) )./(1-pi(j,d) + (r-1).* theta(j,d)).^2;   
        end
       end
       fish1=fish1 + posterior(i,j) .* abs(sum_X1 - sum_Y2);     
       fish2=fish2 + posterior(i,j) .* abs(sum_X3 + sum_Y4 - sum_Y1);  
       fish3=fisher{j,d}(1,2) + posterior(i,j) .* abs(sum_X2 - sum_Y3);  
       
       if (isnan(fish1)==1)          
          break;
      end
     end      
       fisher{j,d}(1,1)= fish1;     
       fisher{j,d}(2,2)= fish2 ;
       fisher{j,d}(1,2)= fish3;  
       fisher{j,d}(2,1)= fisher{j,d}(1,2);  
  end
end

