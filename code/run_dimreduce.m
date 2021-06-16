%% Example 2:  Extract single-trial neural trajectories from Raw Spike Trains using DimReduce
%  Perform dimensionality reduction with DimReduce on raw spike trains 
%  and automatically visualize the extracted single-trial neural 
%  trajectories in DataHigh.  The data has 56 trials for two reach
%  directions. Data is from (Yu et al., 2009).
%  Please see User Guide and online videos for detailed instructions. 
%
%  Quick start (with GPFA trick):
%  1. 20ms bin width
%     1.0 spikes/s  threshold
%     (unchecked) trial-averaged neural trajs
%  2. Method: GPFA
%  3. Select a dimensionality of 40.
%  5. Click the 'Perform Dim Reduction' button (will take ~2min).
%  6. When the PostDimReduce figure pops up, select a dimensionality of 8.
%  7. Click 'Upload to DataHigh.'
%  8. See ex2_visualize.m to visualize the single-
%     trial neural trajectories with DataHigh.
%

epoch = 'f1';
pulse = 16;
angle = 0;

stuff = '/Users/kendranoneman/Projects/mayo';
cd ..

% pulseVnopulse
%load(sprintf('%s/data/eachepoch/pulse-nopulse/rawspiketrains-%s-p%d-d%d.mat',stuff,epoch,pulse,angle));

% fourdirs
%load(sprintf('%s/data/eachepoch/fourdirs/rawspiketrains-%s-p%d.mat',stuff,epoch,pulse));

% pulsespeeds
%load(sprintf('%s/data/eachepoch/pulsespeeds/rawspiketrains-allp-%s-d%d.mat',stuff,epoch,angle));

% manifold (BROKEN)
load(sprintf('%s/data/manifold/datahigh/pall-ef1-d000.mat',stuff));

cd DataHigh
% D(itrial).data : (num_neurons x num_1ms_bins)
DataHigh(D, 'DimReduce');
%cd ./examples