function [accuracy, true_labels, CM, precision,recall,F] = calculateAccuracy(yte, y,classes)
%# Function for calculating clustering accuray and matching found 
%# labels with true labels. Assumes yte and y both are Nx1 vectors with
%# clustering labels. Does not support fuzzy clustering.
%#
%# Algorithm is based on trying out all reorderings of cluster labels, 
%# e.g. if yte = [1 2 2], try [1 2 2] and [2 1 1] so see witch fit 
%# the truth vector the best. Since this approach makes use of perms(),
%# the code will not run for unique(yte) greater than 10, and it will slow
%# down significantly for number of clusters greater than 7.
%#
%# Input:
%#   yte - result from clustering (y-test)
%#   y   - truth vector
%#
%# Output:
%#   accuracy    -   Overall accuracy for entire clustering (OA). For
%#                   overall error, use OE = 1 - OA.
%#   true_labels -   Vector giving the label rearangement witch best 
%#                   match the truth vector (y).
%#   CM          -   Confusion matrix. If unique(yte) = 4, produce a
%#                   4x4 matrix of the number of different errors and  
%#                   correct clusterings done.

N = length(y);
cluster_names = unique(yte);
accuracy = 0;
maxInd = 1;

perm = perms(unique(y));
[pN pM] = size(perm);

true_labels = y;

for i=1:pN
    flipped_labels = zeros(1,N);
    for cl = 1 : pM
        flipped_labels(yte==cluster_names(cl)) = perm(i,cl);
    end

    testAcc = sum(flipped_labels == y')/N;
    if testAcc > accuracy
        accuracy = testAcc;
        maxInd = i;
        true_labels = flipped_labels;
    end

end

CM = zeros(pM,pM);
for rc = 1 : pM
    for cc = 1 : pM
        CM(rc,cc) = sum( ((y'==rc) .* (true_labels==cc)) );
    end
end
confus=CM;
precision=[];
recall=[];
F=[];


for i=1:length(classes)
    S = sum(confus(i,:));
    n=sum(confus(:,i));
    if nargout>=4
        if S
            recall(i) = confus(i,i) / sum(confus(i,:));
        else
            recall(i) = 0;
        end
    end
    S =  sum(confus(:,i));
    if nargout>=3
        if S
            precision(i) = confus(i,i) / S;
         
        else
            precision(i) = 0;
        end
    end
    if nargout>=5
        if (precision(i)+recall(i))
            F(i) = 2 * (precision(i)*recall(i)) / (precision(i)+recall(i));
        else
            F(i) = 0;
        end
    end
    
    
%%compute precision and recall, average precision
% [tp, fp, p] = vl_tpfp(true_labels(1:n), y(1:n)) ;
% 
% % compute precision and recall
% small = 1e-10 ;
% recall = tp / max(p, small) ;
% 
% precision = max(tp, small) ./ max(tp + fp, small) ;
% % 
% % 
% % 
% %   % average precision (for each recalled positive sample)
%    sel = find(diff(recall)) + 1 ;
%    ap(i) = sum(precision(sel)) / p ;
% % 
% 
score=(true_labels(1:n)==y(1:n)');
ap (i)=  sum(precision(i).* (score>0)) / S ;
% % 
 auc (i) = 0.5 * sum((precision(1:end-1) + precision(2:end)) .* diff(recall)) ; 
%   
% [RECALL, PRECISION, INFO] = vl_pr(true_labels(1:n), y(1:n));
% AP(i)=INFO.ap;
end

