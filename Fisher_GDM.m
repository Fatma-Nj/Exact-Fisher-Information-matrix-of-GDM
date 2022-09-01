%% Exact Fisher infromation matrix for Generalized Dirichlet Multinomial distributions

%load data

load EmoBank_bag_of_30_words_kdtree
X=feature_hist;
% Cross varidation (train: 70%, test: 30%)
cv = cvpartition(size(X,1),'HoldOut',0.2);
idx = cv.test;
% Separate to training and test data
dataTrain = X(~idx,:);
dataTest  = X(idx,:);


K=2;  %number of components
Kmax=2; % change it for model slection


%% Fisher scoring estimation algorithm of GDM
tstart = cputime;
[user,sys] = memory;%memory in the begging of your program
[p, pi, theta,MML, MDL,AIC]=fisher_scoring(X,Kmax);
tend = cputime - tstart;
f = @() fisher_scoring(X,Kmax);
timeit(f)

[user2,sys2] = memory;%memory in the end of your program
memory_used_in_bytes=user2.MemAvailableAllArrays-user.MemAvailableAllArrays;


 load Emobank_label


actual_label = emobank_label;
Pc=clusteringGDM(X,K,abs(pi),abs(theta),p);

[pv0,predict_label]=max(Pc,[],2);

%  r= corr(predict_label,actual_label,'Type','Pearson'); % pearson correlation
%  s= corr(predict_label,actual_label,'Type','Spearman') ;% spearman rank correlation
 classes =[1:max(max(actual_label),max(predict_label))];
 
 [accuracy, true_labels, CM, precision,recall,F] = calculateAccuracy(actual_label,predict_label,classes);
  accuracy
  mean(precision)
  mean(recall)
  mean(F)
  

 [a,b]=hist(actual_label,unique(actual_label));
 [c,d]=hist(predict_label,unique(predict_label));






