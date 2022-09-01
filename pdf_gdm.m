function pdf_GDM = pdf_gdm(X,K,pi,theta)
%Calculate the probability density function of the Generalized Dirichlet
%multinmial with the new parametrization of (pi=alpha/alpha+beta ,
%theta=1/alpha+beta)


[N,D]=size(X);
pdf_GDM=size(N,K);
m=sum(X,2);
for i=1:N
   for j=1:K 
     sum_X=0;  sum_Y=0;sum_r=0; sum_tot=0;
     for d=1:D-1
        Y=sum(X(i,d:D));
        Y_l=sum(X(i,d+1:D));
            
       for r=1:X(i,d)
        sum_X = sum_X + log(pi(j,d) + (r-1).* theta(j,d));
       end

       for r=1:Y_l
        sum_Y = sum_X + log(1 - pi(j,d) + (r-1).* theta(j,d));
       end
       
       for r=1:Y
        sum_r= sum_r + log(1 + (r-1).* theta(j,d));
       end
        sum_tot = sum_tot + (sum_X + sum_Y - sum_r);
     end
     pdf_GDM(i,j)= ( sum_tot + gammaln(m(i)+1) - sum(gammaln(X(i,:)+1)) ) ;
       if (isinf(pdf_GDM(i,j))==1)
         pdf_GDM(i,j)=1e10;
       end
   end
end
pdf_GDM = abs(exp (-1e-6 .* pdf_GDM));
%./ sum(exp(-1e-5 * pdf_GDM),1) ;  %%softmax function (beta=1e-5)

end


