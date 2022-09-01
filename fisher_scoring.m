function [p, pi, theta, MML, MDL,AIC]=fisher_scoring(X,Kmax)

[N,D]=size(X);
d=D-1;
%% Fisher scoring estimation algorithm for Generalized Dirichlet multinomial
MML=zeros(Kmax,1);MDL=zeros(Kmax,1);AIC=zeros(Kmax,1);

%for K=1:1:Kmax
K=Kmax;
%%intiatialization of parameters with method of moment
tau=0.9;
[alpha0,beta0,p0] = initial_MoM(X,K);
%%define the new set of parameters
pi0=(alpha0./(alpha0+beta0))  ;
theta0=(1./(alpha0+beta0))  ;
fisher=  exact_fisher_GDM(X,K,pi0,theta0,p0,tau);
Theta_0=cell(K,d);
for j=1:K
    for h=1:d
        Theta_0{j,h}=[pi0(j,h); theta0(j,h)] ;
    end
end


iter_max=2;
p=zeros(1,K);
score_pi=zeros(K,d);
score_theta=zeros(K,d);

posterior=posterior_GDM(X,K,pi0,theta0,p0,tau);
sum_X1=0; sum_X2=0; sum_Y1=0; sum_Y2=0; sum_Y3=0;
loglike = loglikelihood_GDM(X,K,pi0,theta0,p0,tau);

%while(tau<=1) 
    
for iter=1:iter_max
iter




  Score=cell(K,d); Theta=cell(K,d);
    %mixing weight
     p=1/N .* (sum(posterior,1));
     
  for j=1:K      
    for h=1:d
      for i=1:N    
        Y=sum(X(i,h:D));
        Y_l=sum(X(i,h+1:D));
       for r=1:X(i,h)
        sum_X1 = sum_X1 + 1./(pi0(j,h) + (r-1).* theta0(j,h));
        sum_X2 = sum_X2 + (r-1)./(pi0(j,h) + (r-1).* theta0(j,h));
       end  
       
       for r=1:Y_l
        sum_Y1 = sum_Y1 + 1./(1-pi0(j,h) + (r-1).* theta0(j,h));
        sum_Y2 = sum_Y2+ (r-1)./(1-pi0(j,h) + (r-1).* theta0(j,h));        
       end  
       
       for r=1:Y
        sum_Y3 = sum_Y3+ (r-1)./(1 + (r-1).* theta0(j,h));
       end
       %%calculating the score vectors
       score_pi(j,h) = score_pi(j,h) + posterior(i,j) .* abs(sum_X1 - sum_Y1);
       
       score_theta(j,h) = score_theta(j,h) + posterior(i,j) .* abs(sum_X2 + sum_Y2- sum_Y3);
       
     end      
     
       
       
       
      
       Score{j,h}=[score_pi(j,h); score_theta(j,h)];
       %%Fisher update
     
       Theta{j,h}= Theta_0{j,h} + fisher{j,h}\ Score{j,h};
       
       %%updating the parameters
       pi0(j,h) = Theta{j,h}(1,1);    
       theta0(j,h) = Theta{j,h}(2,1);
      
       Theta_0{j,h}= Theta {j,h};
     end             
  end

 % %%smoothing the parameters
 pi0=(pi0 * 1e1)./ ( sum(pi0*1e1,1));
 theta0=(theta0*1e1) ./ ( sum(theta0*1e1,1));
%  
 posterior=posterior_GDM(X,K,pi0,theta0,p,tau);
 fisher=  exact_fisher_GDM(X,K,pi0,theta0,p,tau);
 %  %% normalizing
%  pi0=exp(pi0)./ sum(exp(pi0));
%  theta0=exp(theta0)./sum(exp(theta0));
 %  %%convergence test (maximizing log-likelihood)
   Log_like = loglikelihood_GDM(X,K,pi0,theta0,p,tau);
  if ( (Log_like - loglike) < 1e-3 )
      break;
  end

   
 
end

 tau=tau*5;

 pi=pi0;
 theta=theta0;
%%====================check the MMl criterion============================
%%beta prior for the parameters
alpha1=0.1 .* ones(K,d); beta1=0.02.* ones(K,d);
alpha2= 0.2 .* ones(K,d); beta2=0.04.* ones(K,d);
%beta prior
% prior_theta = gamma(alpha1 + beta1)./ (gamma(alpha1) .* gamma(beta1)) .* (theta.^(alpha1-1)) .* (1-theta).^(beta1-1) ;
% prior_pi = gamma(alpha2 + beta2)./ (gamma(alpha2) .* gamma(beta2)) .* (pi.^(alpha2-1)) .* (1-pi).^(beta2-1) ;
%gamma prior
prior_theta = gamrnd(beta1,alpha1);
%(beta1.^alpha1) .* theta.^(alpha1-1) .* exp(-beta1.*theta)./gamma(alpha1);
prior_pi    = gamrnd(beta2,alpha2);
%(beta2.^alpha2) .* pi.^(alpha2-1) .* exp(-beta2.*pi)./gamma(alpha2);
prod_fisher=1;
for j=1:K
    for h=1:d
        prod_fisher=prod_fisher .* det(fisher{j,h}) ;
    end
end
prod_fisher = prod_fisher .* N./(prod(p)); 
prior = factorial(K-1) .* prod(prod(prior_theta.* prior_pi,1),2) ;
MML (K)= - Log_like - log(prior+1e-10) + 0.5 * log(prod_fisher+1e-10) + 0.5 * K * (2 * D + 1) * ( 1 + log(1/12) );
MDL (K)= - Log_like + 0.5 * K * (2 * D + 1) .* log(N);
AIC (K)= - Log_like + 0.5 * K * (2 * D + 1);


%end

end

%end


