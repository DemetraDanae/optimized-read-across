%% Rulette construct function

function ChromScoreNormAcc = ruletteConstruct(ChromScore)

nChrom = size(ChromScore ,1);

% Calculate the sum of the scores
ScoreSum = sum(ChromScore (:,2));

% Normalize the scores between 0 and 1
ChromScoreNorm = [ ChromScore(:,1) , ChromScore(:,2) ./ScoreSum ];

% Sort ChromScoreNorm rows by 2nd column (score), descending 
%DEACTIVATED (save computational time)
%ChromScoreNorm = sortrows(ChromScoreNorm ,-2);

% Accumulated score
clear ChromScoreNormAcc
for i=1: nChrom
    ChromScoreNormAcc(i,1) = ChromScoreNorm(i,1);
    if (i==1)
        ChromScoreNormAcc(i,2) = ChromScoreNorm(i,2);
    else
        ChromScoreNormAcc(i,2) = ChromScoreNormAcc(i-1,2) + ChromScoreNorm(i,2);
    end
end

end