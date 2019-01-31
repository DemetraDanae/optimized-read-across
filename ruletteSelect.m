%% Rulette select function

function selectedChrom = ruletteSelect(ChromScoreNormAcc)
nChrom = size(ChromScoreNormAcc ,1);

% Generate a random number between 0 and 1
R = rand(1);

% Select the first chromosome with score greater than R
for i=1: nChrom
    if ( ChromScoreNormAcc(i,2) > R )
        selectedChrom = ChromScoreNormAcc(i,1);
        break
    end
end

end