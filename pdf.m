function pdf_GDM = pdf(X,K,pi,theta)
%Calculate the probability density function of the Generalized Dirichlet
%multinmial with the new parametrization of (pi=alpha/alpha+beta ,
%theta=1/alpha+beta)


[N,D]=size(X);
pdf_GDM=size(N,K);
m=sum(X,2);
prod_X=1; prod_Y=1; prod_r=1; prod_tot=1; 
for i=1:N
   for j=1:K       
     for d=1:D
        Y=sum(X(i,h:D));
        Y_l=sum(X(i,h+1:D));
       for r=1:X(i,d)
        prod_X = prod_X .* (pi(j,d) + (r-1).* theta(j,d));
       end
       for r=1:Y_l
        prod_Y = prod_Y .* (1-pi(j,d) + (r-1) .* theta(j,d));
       end
       for r=1:Y
        prod_r= prod_r .* (1 + (r-1).* theta(j,d));
       end
        prod_tot = prod_tot .* ((prod_X .* prod_Y)./prod_r);
      end
   end
       pdf_GDM(i,j)=factorial(m)./prod(factorial(X(i,:))) .* prod_tot;
end
end

