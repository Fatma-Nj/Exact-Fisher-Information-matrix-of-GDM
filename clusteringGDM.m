function [ Pc ] = clusteringGDM (X,K,pi,theta,p)


 pdf=pdf_gdm(X,K,pi,theta);
           
 %pdf_DM=pdf_DM +1e20;         
 Pc=bsxfun(@times,pdf,p);



end

