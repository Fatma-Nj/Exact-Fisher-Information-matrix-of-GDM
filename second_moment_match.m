function [alpha, beta] =second_moment_match(X)
%% estimating the generalized dirichlet parameters using first and second moment

[~,D]=size(X);
d=D-1; 
alpha=zeros(d,1); beta=zeros(d,1);
mu=zeros(d,1); S=zeros(d,1); A=zeros(d,1); B=zeros(d,1); a=zeros(d,1);

mu(1)=mean(X(:,1))+1e-10;
S(1)=var(X(:,1))+1e-10;
A(1)=1; B(1)=1; a(1)=1;
alpha(1)= ( (mu(1).^2) - mu(1) .* (S(1) + mu(1).^2) ) ./ ((S(1) + mu(1).^2) - mu(1).^2) ;
beta(1) = alpha(1) .* ( 1- mu(1) ) ./ mu(1);
for h=2:d 
    mu(h)=mean(X(:,h))+1e-10;
    S(h)=var(X(:,h)) +1e-10;
    

    a(h)=mu(h)./A(h-1);
    
    alpha(h) = ( (a(h) .* mu(h) .* B(h-1)) - mu(h) .* (S(h) + mu(h).^2) ) ./ (A(h-1) .* (S(h) + mu(h).^2) - a(h) .* mu(h).* B(h-1) +1e-10 ) ;
    beta(h) = alpha(h) .* ( A(h-1)- mu(h) ) ./ mu(h) +1e-10;
    prod_a=1; prod_b=1;
    for j=1:h
     prod_a = prod_a .* ( beta(j) ./(alpha(j) + beta(j)));
     prod_b = prod_b .* ( (beta(j) .* ( beta(j)+1))./ ((alpha(j) + beta(j)) .* (alpha(j) + beta(j)+1)));
    end
    A(h) = prod_a;
    B(h) = prod_b;
end    
alpha=abs(alpha); beta=abs(beta);
end