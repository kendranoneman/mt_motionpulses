clear vars;

%% Load data
stuff = '/Users/kendranoneman/Projects/mayo';
loc = 'data/manifold/datahigh'; loc2 = 'data/manifold/counts';
loc3 = 'data/manifold/spikes';
d = load('-mat',sprintf('%s/data/combinedMaestroSpkSortFEF(pa59pulsA).mat',stuff));

folder = sprintf('%s/%s',stuff,loc); folder2 = sprintf('%s/%s',stuff,loc2);
folder3 = sprintf('%s/%s',stuff,loc3);
%% Setup functions
max_neurons = good_neurons(d);

window = 200;
trialvunit(d,window,max_neurons,stuff); 

%% Loading data & making struct S
angles = [0 90 180 270];                                                                                          
epochs = ['f1'; 'f2'; 'p1'; 'p2'];  

for ang = 1:size(angles,2)
    for epo = 1:size(epochs,1)
        S = struct([]);  
        angle = angles(ang);   
        epoch = epochs(epo,:);       
        
        %condition name
        filePattern = fullfile(folder, sprintf('*-e%s-d%3.3d.mat',epoch,angle));
        filePattern2 = fullfile(folder2, sprintf('*-e%s-d%3.3d.mat',epoch,angle));
        filePattern3 = fullfile(folder3, sprintf('*-e%s-d%3.3d.mat',epoch,angle));
        files = dir(filePattern); files2 = dir(filePattern2);  
        files3 = dir(filePattern3);

        for i=1:length(files) 
            S(i).name = files(i).name; 
            S(i).speed = files(i).name(2:3);  
            S(i).epoch = files(i).name(6:7);  
            S(i).angle = files(i).name(10:12);   
            S(i).spikes = load('-mat', sprintf('%s/%s',folder3,files3(i).name));
            S(i).counts = load('-mat', sprintf('%s/%s',folder2,files2(i).name)); 
            S(i).datahigh = load('-mat', sprintf('%s/%s',folder,files(i).name));

        end

        combo = [S(1).counts.N_counts
                 S(2).counts.N_counts
                 S(3).counts.N_counts
                 S(4).counts.N_counts
                 S(5).counts.N_counts]; 
        norm_mu = mean(combo,1); % mean spike count for each unit
        norm_std = std(combo,1); % stddev for each unit

        for ii=1:length(files)
            Z_tot = zeros(size(S(ii).counts.N_counts));
            for a=1:size(S(ii).counts.N_counts,1) % each row (trial)
                for b=1:size(S(ii).counts.N_counts,2) % each col (unit)
                    Z_tot(a,b) = (S(ii).counts.N_counts(a,b) - norm_mu(b))/norm_std(b);
                end
            end
            S(ii).mean = norm_mu;
            S(ii).std = norm_std;
            S(ii).norm = Z_tot; 
        end
        save(sprintf('%s/data/manifold/norm/e%s-d%3.3d.mat',stuff,epoch,angle),'S')
    end
end

%% Combining data into single z-scored "direction" group of trials
N = struct([]);
folder = sprintf('%s/data/manifold/norm',stuff);
epochs = ['f1'; 'f2'; 'p1'; 'p2'];
for epo=1:size(epochs,1)
    epoch = epochs(epo,:);
    filePattern = fullfile(folder, sprintf('*e%s-*.mat',epoch));
    files = dir(filePattern);

    count = 0;
    for i=1:length(files)
        load('-mat', sprintf('%s/%s',folder,files(i).name));
        for j=1:size(S,2)
            count=count+1;
            N(count).name = S(j).name;
            N(count).speed = S(j).speed;
            N(count).epoch = S(j).epoch;
            N(count).angle = S(j).angle;
            N(count).counts = S(j).counts.N_counts;
            N(count).datahigh = S(j).datahigh;
            N(count).norm = S(j).norm;
        end
    end
    save(sprintf('%s/data/manifold/norm/all/e%s.mat',stuff,epoch),'N')
end