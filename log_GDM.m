function [log_lik] = log_GDM(X,pi,theta)
%% this function to calculate the log-likelihood of GDM distributions
[N,D]=size(X);
log_1=0; log_2=0; log_3=0;

for i=1:N
   for h=1:D-1
        Y=sum(X(i,h:D));
        Y_l=sum(X(i,h+1:D));
     for r=1:X(i,h)
       log_1=log_1+log(pi(h)+(r-1).* theta(h));
     end
     for r=1:Y_l
         log_2=log_2+log(1-pi(h)+(r-1).* theta(h));
     end
     for r=1:Y
         log_3=log_3+log(1+(r-1).*theta(h));
     end
     log_lik(i,h)=log_1+log_2-log_3;
   end   
end
log_lik=sum(sum(log_lik,1),2);
