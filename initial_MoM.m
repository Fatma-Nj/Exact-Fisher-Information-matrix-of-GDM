function [alpha,beta,phi] = initial_MoM(X,K)

% X=normalize_norm_mean(X);


%X=normalizing(X);
[L1] = kmeans(X, K);

[N,D] = size(X);
d=D-1;
alpha=zeros(K,d); beta=zeros(K,d);
%Data belong to each cluster
f=(1:K);
F=cell(K,1);
for j=1:K
    F{j}=[];
end
cpt=zeros(K,1);
for i=1:N
  for j=1:K
   h=f(j);
     if (L1(i)==h)
         cpt(j)=cpt(j)+1;
         cc=cpt(j);
         F{h}(cc,:)=X(i,:);  
     end
  end
end

for j=1:K
    Xp=F{j} ;
    [alpha(j,:), beta(j,:)] =second_moment_match(Xp);
 
end

 phi=ones(1,K)/K;



end