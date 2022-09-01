function [Lok_like] = loglikelihood_GDM(X,K,pi,theta,p,tau)

pdf=pdf_gdm(X,K,pi,theta);
posterior=posterior_GDM(X,K,pi,theta,p,tau);

Lok_like=sum(sum(posterior .* log(p .* pdf +1e-10),1),2); 
end
