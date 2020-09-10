
%This program is free software: you can redistribute it and/or modify it under the terms 
%of the GNU General Public License as published by the Free Software Foundation, either 
%version 3 of the License, or (at your option) any later version.

%This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
%without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%You should have received a copy of the GNU General Public License along with this program.  
%If not, see here: http://www.gnu.org/licenses/.

clear all

%% Comment in when using Octave
%pkg load statistics

%% Import files
data_file = 'file_name.csv';
data = readtable(data_file,'Delimiter',';','ReadRowNames',true); % Comment out when using Octave
%data = dlmread(data_file,";",1,1); % Comment in when using Octave

% Attribute values (properties)
Instances = data(1:end,2:end);     

% Last column of the descriptors
col = 40;

% Endpoint experimental values
Labels = data(1:end,1);
Labels = table2array(Labels); %Comment out when using Octave

totalInst = size(Instances ,1);
nAttr = size(Instances ,2);

%% Scaling 0-1
Instances = table2array(Instances);
for j=1:nAttr
   for i=1:totalInst
      Instances_sc(i,j) = (Instances(i,j)-min(Instances(:,j)))/(max(Instances(:,j))-min(Instances(:,j)));
   end
end

% Check if descriptors have NA values
flag = 0;
for i=1:col
  if isnan(Instances_sc(:,i))== true
    flag = flag+1;
  end
end
col = col-flag; %New column according to reduced set

Instances_sc = Instances_sc(:,all(~isnan(Instances_sc)));

nAttr = size(Instances_sc ,2);

%% Partitioning

data_sc = [Labels, Instances_sc]; %Initial data normalized

% LOO cross-validation
%per = round(0.66*totalInst);
train_data = data_sc; %(kenstone(data_sc(:,2:end),per),:); %kennard-stones

%test_data = setdiff(data_sc,train_data,'rows');

% Train instances
train = train_data(:,2:end);
nInst = size(train,1);

% Train labels
trainL = train_data(:,1);

% Test instances
%test = test_data(:,2:end);
%nInst_test = size(test,1);

% Test labels
%testL = test_data(:,1);

%% Parameters

% Number of chromosomes for each generation (an even number)
nChrom = 100;

% Probability for a gene to have value "1" intitially
initGeneProb = 0.60;

% Probability for mutation of each gene
mutProb = 0.01;

% Probability for crossover between chromosomes
crossProb = 0.70;

% Non uniform mutation probability
nonUnf = 0.1;

for i=1:nInst
  for j=1:nInst
    Dist(i,j) = norm(train(i,:)-train(j,:));  %Euclidean
  end
end

% Lower bound of the threshold
lGA_min = 0.1;

% Upper bound of the threshold
lGA_max = mean(max(Dist));

% Freezing parameter
bSA = 1;

% Maximum number of generations after the initial
maxGenerations = 1000;

% Number of samples with a prediction
satSamples = round(0.3*nInst);

%% Initial population 

% Create a random initial population of chromosomes
Chromosomes = [round(rand(nChrom ,nAttr)), datasample(datasample(Dist,1),nChrom)']; %Attributes-Threshold 1 comment out when using Octave
%Chromosomes = [round(rand(nChrom ,nAttr)), rand(1,nChrom)'];  %Comment in when using Octave

for i=1:nChrom
    for j=1:nAttr
        a = rand(1);
        if ( a>initGeneProb )
            Chromosomes(i,j) = 0;
        else
            Chromosomes(i,j) = 1;
        end
    end
end

% Test the corresponding models
ch = 1; 
while ch <= nChrom
    for j=1: nInst
        LocalInstances(j,:) = train(j,:) .*Chromosomes(ch,1:nAttr);
    end
           
    ChromScore(ch,1) = ch;
    [ChromScore(ch,2), predicted] = read(LocalInstances, trainL, Chromosomes(ch,nAttr+1));
    
    if ((satSamples > predicted)  || (ChromScore(ch,2) < 0.001))
        Chromosomes(ch,:) = [round(rand(1 ,nAttr)), datasample(datasample(Dist,1),1)]; %comment out when using Octave
        %Chromosomes(ch,:) = [round(rand(1 ,nAttr)), rand(1,1)']; %comment in when using Octave
        for j=1:nAttr
          a = rand(1);
          if ( a>initGeneProb )
            Chromosomes(ch,j) = 0;
          else
            Chromosomes(ch,j) = 1;
          end
        end
        
        % Also check if there is at least one selected variable
        if sum(Chromosomes(ch,1:nAttr))==0
            randomGene = datasample(1:nAttr,1); %Make one random gene equal to 1
            Chromosomes(ch,randomGene) = 1;
        end
        
    else
        ch = ch+1;
    end
    clear LocalInstances; 
end

% Best score
max_score = ChromScore(1,2);
bestChrom = Chromosomes(1,:);
bestGen = 1;
for i=2:nChrom
    if (ChromScore(i,2)>max_score)
        max_score = ChromScore(i,2);
        bestChrom = Chromosomes(i,:); %Selected variables and thresholds
    end
end

'Initial'

%% Main loop

for generation =1:maxGenerations
     
    % Create the rulette for selection (one time)
    % Accumulated Normalized Score for Chromosomes - an (nChrom ,2) array
    ChromScoreNormAcc = ruletteConstruct(ChromScore);

    counter = 0;
    while counter <= nChrom -2
    % counter: how many new chromosomes have been produced?

        % Select two chromosomes
        i1 = ruletteSelect(ChromScoreNormAcc);
        i2 = ruletteSelect(ChromScoreNormAcc);
        
        % Best chromosome between the 2 selected ones
        bestParent = i1;
        if (ChromScore(i1,2)<ChromScore(i2,2))
            bestParent = i2;
        end
        
        % Crossover probability
        cr = rand(1);
        if ( cr<=crossProb )
        % The chromosomes will be crossovered
            
            % Set a random point for crossover
            point = floor(rand(1)*nAttr);
            % Reproduce the selected chromosomes
            newChrom = reproduce(point,i1,i2,Chromosomes);
         else
            newChrom(1,:) = Chromosomes(i1,:);
            newChrom(2,:) = Chromosomes(i2,:);
        end

        % Mutations
        for i=1:2
            for j=1:nAttr
            % Toggle the value of a gene with a probability mutProb
                if ( rand(1) < mutProb )
                    if (newChrom(i,j) == 0)
                        newChrom(i,j) = 1;
                    else
                        newChrom(i,j) = 0;
                    end % end if newChrom
                 end % end if rand(1)
             end % end for j
             
             % Also check if there is at least one selected variable
            if sum(newChrom(i,1:nAttr))==0
               randomGene = datasample(1:nAttr,1); %Make one random gene equal to 1
               newChrom(i,randomGene) = 1;
            end
             
             % Mutation of threshold
             if ( rand(1)< nonUnf )
                r = rand(1); 
                if randi([0 1])==0
                    newChrom(i,nAttr+1) = newChrom(i,nAttr+1)+(lGA_max-newChrom(i,nAttr+1))*(1-r^(1-generation/maxGenerations)^bSA); 
                else
                    newChrom(i,nAttr+1) = newChrom(i,nAttr+1)-(newChrom(i,nAttr+1)-lGA_min)*(1-r^(1-generation/maxGenerations)^bSA); 
                end
             end
        end % end for i

        % Test the new chromosomes and store the scores
        for i=1:2
            counter = counter + 1;
            newChromScore(counter ,1) = counter;
            for j=1:nInst
                LocalInstances(j,:) = train(j,:).*newChrom(i,1:nAttr);
            end
            
            [newChromScore(counter ,2), predicted] = read(LocalInstances, trainL, newChrom(i,nAttr+1));
            
            if (satSamples-predicted)>0 
                newChrom(i,:) = Chromosomes(bestParent,:);
                newChromScore(counter,2) = ChromScore(bestParent,2);
            end
            
            clear LocalInstances; 
        end
                        
        % Add new chromosomes to the new generation
        if ( counter == 2 )
            newGeneration = newChrom;
        else
            newGeneration = [ newGeneration ;newChrom ];
        end

        % Clear the (temporary) newChrom
        clear newChrom
    end

    % Replace the old generation with the new
    Chromosomes = newGeneration;
    ChromScore = newChromScore;
    
    % Insert best chromosome into the population (if not already included)
    % elitism 
    if any(ChromScore(:,2) > max_score) == 0
       pos = find(ChromScore(:,2) == min(ChromScore(:,2)));
       ChromScore(pos(1),2) = max_score;
       Chromosomes(pos(1),:) = bestChrom;
    end 
    
    % Best score
    for ch=1:nChrom
        if (ChromScore(ch,2) > max_score)
            max_score = ChromScore(ch,2);
            bestChrom = Chromosomes(ch,:);
            bestGen = generation;
        end
    end

    % Clear the new generation
    clear newGeneration
    clear newChromScore
    
    generation
end

genome = bestChrom;

% Test the dominant model
for j=1:nInst
    LocalInstances(j,:) = train(j,:) .* genome(1:nAttr);
end

%for j=1:nInst_test
%    LocalInstances_test(j,:) = test(j,:) .* genome(1:nAttr);
%end

% Distances & neighbours
clear Dist
for i=1:nInst
    for j=1:nInst
        Dist(i,j) = norm(LocalInstances(i,:)-LocalInstances(j,:));  %Euclidean
        if (Dist(i,j)<=genome(nAttr+1)) %Neighbours selection
            neib(i,j) = 1;
        else
            neib(i,j) = 0;
        end
        neib(i,i) = 0;
    end
end

% Weighting factors
wf = 1./(1+Dist);
       
RAV = ((neib.*wf)*trainL)./sum(neib.*wf,2); %LOO CV
     
final = [trainL RAV];
final = final(all(!isnan(final),2),:);
       
MSE = 1/size(final,1)*sum((final(:,1)-final(:,2)).^2);
     
scatter(final(:,1),final(:,2))

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

%% Tropsha's acceptability criteria
[verdict, R2] = tropsha(final)

thr = genome(nAttr+1)
var = sum(genome(1:nAttr))

%% Results final
res(1,1) = R2; %R2
res(1,2) = 1-sum((final(:,1)-final(:,2)).^2)/sum((final(:,1)-mean(trainL)).^2); %external%q2
res(1,3) = 1/MSE; %score
res(1,4) = size(final,1); %predicted samples
res(1,5) = bestGen;
res(1,6) = thr;
res(1,7) = var;
