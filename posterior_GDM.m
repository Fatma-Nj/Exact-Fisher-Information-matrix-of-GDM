function posterior=posterior_GDM(X,K,pi,theta,p,tau)


pdf_GDM=pdf_gdm(X,K,pi,theta);
  
posterior=(p .* pdf_GDM).^tau./ (sum (p.* pdf_GDM,2).^tau);

end

