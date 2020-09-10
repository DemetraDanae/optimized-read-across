%% RAV function

function [score, finalSamples] = read(LocalInstances, Labels, Thresholds)
  pkg load statistics
 nInst = size(LocalInstances,1);
    
% Distances
for i=1:nInst
    for j=1:nInst
        Dist(i,j) = norm(LocalInstances(i,:)-LocalInstances(j,:));  %Euclidean
    end
end
    
% Weighting factors
wf = 1./(1+Dist);
    
% Neighbour selection
for ref=1:nInst
    for rem=1:nInst
        if (Dist(ref,rem) < Thresholds)
            neib(ref,rem) = 1;
        else
            neib(ref,rem) = 0;
        end 
    end
    neib(ref,ref) = 0;
end          

RAV = ((neib.*wf)*Labels)./(sum(neib.*wf))'; %LOO CV

final = [Labels RAV];
final = final(all(!isnan(final),2),:);
       
MSE = 1/size(final,1)*sum((final(:,1)-final(:,2)).^2);
    
if MSE==NaN
    error = 1E+10;
else
    error = MSE;
end

if error==1E+10 %no neighbours selected
    score = 0;
else
    OF = error;
    score = 1/(OF+10E-5);
end
    
finalSamples = size(final,1);
end
