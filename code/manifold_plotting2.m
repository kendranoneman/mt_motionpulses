clear all;

stuff = '/Users/kendranoneman/Projects/mayo'; loc = '/data/manifold/norm/all';
green = [51 102 0]/255;
yellow = [245 246 81]/255;

%% Combining directions, plotting points/plane
% X = trials x units (nxp)
% coeff = principal component coefficients (pxp), each col is coefficients
% for one principal component (descending order of component variance)
% score = reps of X in principal component space (trials x components)
% latent = principal component variances (eigenvalues of covariance matrix
% of X)
% tsquared = Hotelling's t-squared statistic for each trial in X
% explained = percentage of total variance explained by each principal comp
% mu = estimated mean of each variable in X

%epochs = ['f1'; 'f2'; 'p1'; 'p2'];
epochs = ['f1'];
mean_error_n123 = []; mean_error_p123 = [];
for j=1:size(epochs,1)
    epoch = epochs(j,:);
    load(sprintf('%s/%s/e%s.mat',stuff,loc,epoch));

    % Splitting into pulse speeds
    d_01 = []; d_02 = [];
    d_04 = []; d_08 = [];
    d_16 = [];
    for i = 1:size(N,2)
       if isequal(N(i).speed,'01')
           d_01 = [d_01; N(i).norm];
       elseif isequal(N(i).speed,'02')
           d_02 = [d_02; N(i).norm];
       elseif isequal(N(i).speed,'04')
           d_04 = [d_04; N(i).norm];
       elseif isequal(N(i).speed,'08')
           d_08 = [d_08; N(i).norm];
       elseif isequal(N(i).speed,'16')
           d_16 = [d_16; N(i).norm];
       end
    end

    speeds = [1, 2, 4, 8, 16];
    for k=1:size(speeds,2)
        speed = speeds(k);

        if speed == 1
            neuron1 = d_01(:,1); neuron2 = d_01(:,2); neuron3 = d_01(:,3);
            num_trials = size(d_01,1);
            [coeff,score,latent,tsquared,explained,mu] = pca(d_01);
        elseif speed == 2
            neuron1 = d_02(:,1); neuron2 = d_02(:,2); neuron3 = d_02(:,3);
            num_trials = size(d_02,1);
            [coeff,score,latent,tsquared,explained,mu] = pca(d_02);
        elseif speed == 4
            neuron1 = d_04(:,1); neuron2 = d_04(:,2); neuron3 = d_04(:,3);
            num_trials = size(d_04,1);
            [coeff,score,latent,tsquared,explained,mu] = pca(d_04);
        elseif speed == 8
            neuron1 = d_08(:,1); neuron2 = d_08(:,2); neuron3 = d_08(:,3);
            num_trials = size(d_08,1);
            [coeff,score,latent,tsquared,explained,mu] = pca(d_08);
        elseif speed == 16
            neuron1 = d_16(:,1); neuron2 = d_16(:,2); neuron3 = d_16(:,3);
            num_trials = size(d_16,1);
            [coeff,score,latent,tsquared,explained,mu] = pca(d_16);
        end

        [X, Y] = meshgrid(linspace(min(neuron3),max(neuron3)), linspace(min(neuron2),max(neuron2)));
        XYZ_1 = [X(:) Y(:) ];
        plot3(XYZ_1(:,1),XYZ_1(:,2),XYZ_1(:,3),'r.');
        

        % PCA 
%         p1 = score(:,1); p2 = score(:,2); p3 = score(:,3);
%         num_trials2 = size(score,1);
%         A = [p1 p2 ones(num_trials2,1)];
%         b = p3;
%         
%         [X2 Y2] = meshgrid(linspace(min(p1),max(p1)), linspace(min(p2),max(p2)));

    end
end
