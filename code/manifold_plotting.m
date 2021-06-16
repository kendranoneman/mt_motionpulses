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

epochs = ['f1'; 'f2'; 'p1'; 'p2'];
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

    XYZ1 = struct([]);
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

        % Construct matrices (Ax = b)
        A = [neuron1 neuron2 ones(num_trials,1)]; % matrix A
        b = neuron3; % matrix b 

        % Solve in the least squares sense
        fit = inv((A'*A)) * A'*b; % fitting params
        errors = b - A*fit;
        mean_error = mean(errors); % close to zero since some +/-
        mdl = norm(errors,1); % Euclidean norm of errors, sum of element magnitudes
        mean_error_n123 = [mean_error_n123; mdl];
        
%         h1 = figure('Visible','off'); 
%         histogram(errors,'BinWidth',0.5);
%         xlabel('Error (b - A*fit)');
%         ylabel('Number of Trials');
%         xline(0,'k--');
%         xlim([-6 6]);
%         ylim([0 60]);
%         title(sprintf('Goodness of Fit for Each Trial to the Least Squares Plane: Speed = %2.2d deg/s, Epoch = %s',speed,epoch));
%         saveas(h1,sprintf('%s/figs/manifold/neurons123/norm/hist-p%2.2d-e%s.png',stuff,speed,epoch));
        
        %fprintf('Solution: %f x + %f y + %f = 1z\n', fit(1),fit(2),fit(3));

        [X1 Y1] = meshgrid(linspace(min(neuron1),max(neuron1)), linspace(min(neuron2),max(neuron2)));
        Z1 = zeros(size(X1));
        for r=1:size(X1,1)
            for c = 1:size(X1,2)
                Z1(r,c) = fit(1) * X1(r,c) + fit(2) * Y1(r,c) + fit(3);
            end
        end
        
%         syms x y z
%         a = fit(1); b = fit(2); c = 1; d = fit(3);
%         f = a*x + b*y + c*z - d;
%         g1 = gradient(f, [x, y, z]);
      
        [Nx,Ny,Nz] = surfnorm(X1, Y1, Z1);
        nx = Nx(1,1); ny = Ny(1,1); nz = Nz(1,1);
        
        v1 = [nx, ny, nz];
        v1 = v1/norm(v1);

%         v1 = [fit(1), fit(2), fit(3)];
%         v1 = v1/norm(v1);
        
        XYZ(k).speed = speed;
        XYZ(k).X1 = X1;
        XYZ(k).Y1 = Y1;
        XYZ(k).Z1 = Z1;
        XYZ(k).v1 = v1;
        
        % Plotting points & plane for 3 neurons
        f = figure('Visible','off');
        mesh(X1,Y1,Z1,'EdgeColor','none','FaceAlpha','0.5','FaceColor',yellow);
        hold on;
        grid off;
        plot3(neuron1,neuron2,neuron3,'.','Color',green,'MarkerSize',15)
        quiver3(median(median(X1)),median(median(Y1)),median(median(Z1)),v1(1),v1(2),v1(3),'r','LineWidth',3);
        title(sprintf('Pulse Speed = %2.2d deg/s, Epoch = %s',speed,epoch));
        xlabel('Unit 1 Normalized FR');
        ylabel('Unit 2 Normalized FR');
        zlabel('Unit 3 Normalized FR');
%         xlim([min(neuron1)-1 max(neuron1)+1]);
%         ylim([min(neuron2)-1 max(neuron2)+1]);
%         zlim([min(neuron3)-1 max(neuron3)+1]);
%         xlim([-5 10]);
%         ylim([-4 6]);
%         zlim([-4 4]);
        str = {[sprintf('Solution: %f x + %f y + %f = z', fit(1),fit(2),fit(3))],[sprintf('Norm Error: %f',mdl)]};
        dim = [.2 .5 .3 .3];
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        saveas(f,sprintf('%s/figs/manifold/neurons123/norm/p%2.2d-e%s.png',stuff,speed,epoch))
        
        % Percentage of total variance explained by each principal comp
        hold off;
%         f2 = figure('Visible','off');
%         h = bar(explained);
%         pe3 = explained(1)+explained(2)+explained(3);
%         dim = [.5 .5 .3 .3]; str = sprintf('First three components explain: %0.2f%%',pe3);
%         annotation('textbox', dim, 'Str', str,'FitBoxToText','on');
%         xlabel('Principle Component');
%         ylabel('Percentage of Total Variance Explained');
%         title(sprintf('Pulse Speed = %2.2d deg/s, Epoch = %s',speed,epoch));
%         saveas(f2,sprintf('%s/figs/manifold/pca/componentvar/norm/p%2.2d-e%s-dall.png',stuff,speed,epoch));
        
        p1 = score(:,1); p2 = score(:,2); p3 = score(:,3);
        num_trials2 = size(score,1);
        A = [p1 p2 ones(num_trials2,1)];
        b = p3;

        fit = inv((A'*A)) * A'*b;
        errors = b - A*fit;
        mean_error = mean(errors); % close to zero since some +/-
        mdl = norm(errors,1); % Euclidean norm of errors 
        mean_error_p123 = [mean_error_p123; mdl];
        
%         h2 = figure('Visible','off'); 
%         histogram(errors,'BinWidth',0.5);
%         xlabel('Error (b - A*fit)');
%         ylabel('Number of Trials');
%         xline(0,'k--'); 
%         xlim([-6 6]); 
%         ylim([0 60]);
%         title(sprintf('Goodness of Fit for Each Trial to the Least Squares Plane: Speed = %2.2d deg/s, Epoch = %s',speed,epoch));
%         saveas(h2,sprintf('%s/figs/manifold/pca/princ123/norm/hist-p%2.2d-e%s.png',stuff,speed,epoch));
        
        %display(sprintf('Solution: %f x + %f y + %f = z', fit(1),fit(2),fit(3)));

        [X2 Y2] = meshgrid(linspace(min(p1),max(p1)), linspace(min(p2),max(p2)));
        %midpt = 
        Z2 = zeros(size(X2));
        for r=1:size(X2,1)
            for c = 1:size(X2,2)
                Z2(r,c) = fit(1) * X2(r,c) + fit(2) * Y2(r,c) + fit(3);
            end
        end

%         syms x y z
%         a = fit(1); b = fit(2); c = 1; d = fit(3);
%         f = a*x + b*y + c*z - d;
%         g2 = gradient(f, [x, y, z]);
        
        
        [Nx,Ny,Nz] = surfnorm(X2,Y2,Z2);
        nx = Nx(1,1); ny = Ny(1,1); nz = Nz(1,1);
        
        v2 = [nx, ny, nz];
        v2 = (v2/norm(v2));

%         a = fit(1); b = fit(2); c = 1; d = fit(3);
%         v2 = [a, b, c];
%         v2 = v2/norm(v2);
        
        XYZ(k).X2 = X2;
        XYZ(k).Y2 = Y2;
        XYZ(k).Z2 = Z2;
        XYZ(k).v2 = v2; 
        
        save(sprintf('%s/figs/manifold/stackedspeeds/xyz/xyz-e%s.mat',stuff,epoch), 'XYZ');
        
        % Points and plane for 3 principal components
        f3 = figure('Visible','off');
        mesh(X2,Y2,Z2,'EdgeColor','none','FaceAlpha','0.5','FaceColor',yellow);
        hold on;
        grid off;
        plot3(p1,p2,p3,'.','Color',green,'MarkerSize',15)
        %plot3(mean(mean(X2)),mean(mean(Y2)),mean(mean(Z2)),'ro')
        quiver3(mean([min(p1),max(p1)]),mean([min(p2),max(p2)]),mean(mean(Z2)),v2(1),v2(2),v2(3),'r','LineWidth',3,'ShowArrowHead','on');
        xlabel('1st Principal Component (Normalized)');
        ylabel('2nd Principal Component (Normalized)');
        zlabel('3rd Principal Component (Normalized)');
%         xlim([min(score(:,1))-1 max(score(:,1))+1]);
%         ylim([min(score(:,2))-1 max(score(:,2))+1]);
%         zlim([min(score(:,3))-1 max(score(:,3))+1]);
%         xlim([-5 10]);
%         ylim([-4 6]);
%         zlim([-4 4]);
        str = {[sprintf('Solution: %f x + %f y + %f = z', fit(1),fit(2),fit(3))],[sprintf('Norm Error: %f',mdl)]};
        dim = [.2 .5 .3 .3];
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        title(sprintf('Pulse Speed = %2.2d deg/s, Epoch = %s',speed,epoch));
        saveas(f3,sprintf('%s/figs/manifold/pca/princ123/norm/p%2.2d-e%s.png',stuff,speed,epoch));
    end
    n123_avg = mean(mean_error_n123);
    p123_avg = mean(mean_error_p123);

    %% Plotting all 5 planes together
    f4 = figure %('Visible','off');
    m1 = mesh(XYZ(1).X1,XYZ(1).Y1,XYZ(1).Z1,'EdgeColor','none','FaceAlpha','0.5','FaceColor',[76 0 153]/255);
    hold on;
    q1 = quiver3(median(median(XYZ(1).X1)),median(median(XYZ(1).Y1)),median(median(XYZ(1).Z1)),XYZ(1).v1(1),XYZ(1).v1(2),XYZ(1).v1(3),'Color',[76 0 153]/255,'LineWidth',3);

    m2 = mesh(XYZ(2).X1,XYZ(2).Y1,XYZ(2).Z1,'EdgeColor','none','FaceAlpha','0.5','FaceColor',[127 0 255]/255);
    q2 = quiver3(median(median(XYZ(2).X1)),median(median(XYZ(2).Y1)),median(median(XYZ(2).Z1)),XYZ(2).v1(1),XYZ(2).v1(2),XYZ(2).v1(3),'Color',[127 0 255]/255,'LineWidth',3);

    m3 = mesh(XYZ(3).X1,XYZ(3).Y1,XYZ(3).Z1,'EdgeColor','none','FaceAlpha','0.5','FaceColor',[153 51 255]/255);
    q3 = quiver3(median(median(XYZ(3).X1)),median(median(XYZ(3).Y1)),median(median(XYZ(3).Z1)),XYZ(3).v1(1),XYZ(3).v1(2),XYZ(3).v1(3),'Color',[153 51 255]/255,'LineWidth',3);

    m4 = mesh(XYZ(4).X1,XYZ(4).Y1,XYZ(4).Z1,'EdgeColor','none','FaceAlpha','0.5','FaceColor',[178 102 255]/255);
    q4 = quiver3(median(median(XYZ(4).X1)),median(median(XYZ(4).Y1)),median(median(XYZ(4).Z1)),XYZ(4).v1(1),XYZ(4).v1(2),XYZ(4).v1(3),'Color',[178 102 255]/255,'LineWidth',3);

    m5 = mesh(XYZ(5).X1,XYZ(5).Y1,XYZ(5).Z1,'EdgeColor','none','FaceAlpha','0.5','FaceColor',[204 153 255]/255);
    q5 = quiver3(median(median(XYZ(5).X1)),median(median(XYZ(5).Y1)),median(median(XYZ(5).Z1)),XYZ(5).v1(1),XYZ(5).v1(2),XYZ(5).v1(3),'Color',[204 153 255]/255,'LineWidth',3);

    lgd = legend('1','','2','','4','','8','','16','');
    %lgd = legend('1','2','4','8','16');
    lgd.FontSize = 14;
    title(sprintf('Planes for Each Pulse Speed (Epoch = %s)',epoch));
    xlabel('Unit 1 Normalized FR');
    ylabel('Unit 2 Normalized FR');
    zlabel('Unit 3 Normalized FR');
    saveas(f4,sprintf('%s/figs/manifold/stackedspeeds/neurons123/vect/e%s.png',stuff,epoch));


    f5 = figure %('Visible','off');
    mesh(XYZ(1).X2,XYZ(1).Y2,XYZ(1).Z2,'EdgeColor','none','FaceAlpha','0.5','FaceColor',[76 0 153]/255);
    hold on;
    q1 = quiver3(median(median(XYZ(1).X2)),median(median(XYZ(1).Y2)),median(median(XYZ(1).Z2)),XYZ(1).v2(1),XYZ(1).v2(2),XYZ(1).v2(3),'Color',[76 0 153]/255,'LineWidth',3);

    mesh(XYZ(2).X2,XYZ(2).Y2,XYZ(2).Z2,'EdgeColor','none','FaceAlpha','0.5','FaceColor',[127 0 255]/255);
    q2 = quiver3(median(median(XYZ(2).X2)),median(median(XYZ(2).Y2)),median(median(XYZ(2).Z2)),XYZ(2).v2(1),XYZ(2).v2(2),XYZ(2).v2(3),'Color',[127 0 255]/255,'LineWidth',3);

    mesh(XYZ(3).X2,XYZ(3).Y2,XYZ(3).Z2,'EdgeColor','none','FaceAlpha','0.5','FaceColor',[153 51 255]/255);
    q3 = quiver3(median(median(XYZ(3).X2)),median(median(XYZ(3).Y2)),median(median(XYZ(3).Z2)),XYZ(3).v2(1),XYZ(3).v2(2),XYZ(3).v2(3),'Color',[153 51 255]/255,'LineWidth',3);

    mesh(XYZ(4).X2,XYZ(4).Y2,XYZ(4).Z2,'EdgeColor','none','FaceAlpha','0.5','FaceColor',[178 102 255]/255);
    q4 = quiver3(median(median(XYZ(4).X2)),median(median(XYZ(4).Y2)),median(median(XYZ(4).Z2)),XYZ(4).v2(1),XYZ(4).v2(2),XYZ(4).v2(3),'Color',[178 102 255]/255,'LineWidth',3);

    mesh(XYZ(5).X2,XYZ(5).Y2,XYZ(5).Z2,'EdgeColor','none','FaceAlpha','0.5','FaceColor',[204 153 255]/255);
    q5 = quiver3(median(median(XYZ(5).X2)),median(median(XYZ(5).Y2)),median(median(XYZ(5).Z2)),XYZ(5).v2(1),XYZ(5).v2(2),XYZ(5).v2(3),'Color',[204 153 255]/255,'LineWidth',3)

    %lgd = legend('1','2','4','8','16');
    lgd = legend('1','','2','','4','','8','','16','');
    lgd.FontSize = 14;
    title(sprintf('Planes for Each Pulse Speed (Epoch = %s)',epoch));
    xlabel('1st Principal Component (Normalized)');
    ylabel('2nd Principal Component (Normalized)');
    zlabel('3rd Principal Component (Normalized)');
    zlim([-1.5e-15 1.5e-15]);
    saveas(f5,sprintf('%s/figs/manifold/stackedspeeds/pca123/vect/e%s.png',stuff,epoch));

    %% Plotting all 5 normal vectors together
    % 3 neurons
    f6 = figure('Visible','off');

    q1 = quiver3(median(median(XYZ(1).X1)),median(median(XYZ(1).Y1)),median(median(XYZ(1).Z1)),XYZ(1).v1(1),XYZ(1).v1(2),XYZ(1).v1(3),'Color',[76 0 153]/255,'LineWidth',3);
    hold on;
    q2 = quiver3(median(median(XYZ(2).X1)),median(median(XYZ(2).Y1)),median(median(XYZ(2).Z1)),XYZ(2).v1(1),XYZ(2).v1(2),XYZ(2).v1(3),'Color',[127 0 255]/255,'LineWidth',3);
    q3 = quiver3(median(median(XYZ(3).X1)),median(median(XYZ(3).Y1)),median(median(XYZ(3).Z1)),XYZ(3).v1(1),XYZ(3).v1(2),XYZ(3).v1(3),'Color',[153 51 255]/255,'LineWidth',3);
    q4 = quiver3(median(median(XYZ(4).X1)),median(median(XYZ(4).Y1)),median(median(XYZ(4).Z1)),XYZ(4).v1(1),XYZ(4).v1(2),XYZ(4).v1(3),'Color',[178 102 255]/255,'LineWidth',3);
    q5 = quiver3(median(median(XYZ(5).X1)),median(median(XYZ(5).Y1)),median(median(XYZ(5).Z1)),XYZ(5).v1(1),XYZ(5).v1(2),XYZ(5).v1(3),'Color',[204 153 255]/255,'LineWidth',3);

    lgd = legend('1','2','4','8','16');
    lgd.FontSize = 14;
    title(sprintf('Planes for Each Pulse Speed (Epoch = %s)',epoch));
    xlabel('Unit 1 Normalized FR');
    ylabel('Unit 2 Normalized FR');
    zlabel('Unit 3 Normalized FR');
    saveas(f6,sprintf('%s/figs/manifold/stackedspeeds/neurons123/vect/vect-e%s.png',stuff,epoch));

    % 3 principal components
    f7 = figure('Visible','off');
    q1 = quiver3(median(median(XYZ(1).X2)),median(median(XYZ(1).Y2)),median(median(XYZ(1).Z2)),XYZ(1).v2(1),XYZ(1).v2(2),XYZ(1).v2(3),'Color',[76 0 153]/255,'LineWidth',3);
    hold on;
    q2 = quiver3(median(median(XYZ(2).X2)),median(median(XYZ(2).Y2)),median(median(XYZ(2).Z2)),XYZ(2).v2(1),XYZ(2).v2(2),XYZ(2).v2(3),'Color',[127 0 255]/255,'LineWidth',3);
    q3 = quiver3(median(median(XYZ(3).X2)),median(median(XYZ(3).Y2)),median(median(XYZ(3).Z2)),XYZ(3).v2(1),XYZ(3).v2(2),XYZ(3).v2(3),'Color',[153 51 255]/255,'LineWidth',3);
    q4 = quiver3(median(median(XYZ(4).X2)),median(median(XYZ(4).Y2)),median(median(XYZ(4).Z2)),XYZ(4).v2(1),XYZ(4).v2(2),XYZ(4).v2(3),'Color',[178 102 255]/255,'LineWidth',3);
    q5 = quiver3(median(median(XYZ(5).X2)),median(median(XYZ(5).Y2)),median(median(XYZ(5).Z2)),XYZ(5).v2(1),XYZ(5).v2(2),XYZ(5).v2(3),'Color',[204 153 255]/255,'LineWidth',3);

    lgd = legend('1','2','4','8','16');
    lgd.FontSize = 14;
    title(sprintf('Planes for Each Pulse Speed (Epoch = %s)',epoch));
    xlabel('1st Principal Component (Normalized)');
    ylabel('2nd Principal Component (Normalized)');
    zlabel('3rd Principal Component (Normalized)');
    %zlim([-1.5e-15 1.5e-15]);
    saveas(f7,sprintf('%s/figs/manifold/stackedspeeds/pca123/vect/vect-e%s.png',stuff,epoch));
end

