 %% Reproduction function
function newChrom = reproduce(point , i1, i2, Chromosomes)

nPoints = size(Chromosomes ,2);

newChrom (1,:) = [ Chromosomes(i1 ,1:point) Chromosomes(i2,point+1: nPoints) ];
newChrom (2,:) = [ Chromosomes(i2 ,1:point) Chromosomes(i1,point+1: nPoints) ];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Developed by Gerasimos Chourdakis - 2013, 2014 %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
