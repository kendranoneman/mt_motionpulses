clear all;

stuff = '/Users/kendranoneman/Projects/mayo'; loc = '/data/manifold/spikes';
%% Compare normalized to non-normalized counts 
epochs = ['f1'; 'f2'; 'p1'; 'p2'];
angles = [0 90 180 270];
speeds = [1 2 4 8 16];

for epo=1:size(epochs,1)
    for ang=1:size(angles,2)
        for spe=1:size(speeds,2)
            epoch = epochs(epo,:);
            angle = angles(ang);
            speed = speeds(spe);
            
            load(sprintf('%s/%s/p%2.2d-e%s-d%3.3d.mat',stuff,loc,speed,epoch,angle));

            P = struct([]);
            for b=1:49
                nn = [];
                for a=1:size(G,2)
                    nn = [nn; G(a).bin20(b,:)];
                end
                P(b).data = nn;
            end    

            %% Plotting PSTHs for each unit
            mu_norm_tot = [];
            for i=1:size(P,2)
                mu = mean(P(i).data,1)*100; % spikes/s

                % PSTH for each unit
%                 f1 = figure('Visible','off');
%                 x = 10:10:200;
%                 y = mu;
%                 bar(x,y)
%                 hold on;
%                 xlabel('Time (10ms bins)')
%                 ylabel('Spikes/s')
%                 title(sprintf('PSTH for each Unit (%s): Epoch = %s, Pulse Speed = %d deg/s, Angle = %d ',i,epoch,speed,angle))
%                 saveas(f1,sprintf('%s/figs/psth/eachunit/p%2.2d-e%s-d%3.3d-%2.2d.png',stuff,speed,epoch,angle,i));

                % Normalize by maximum firing rate
                max_fire = max(mu); % max fr for each unit
                mu_norm_tot = [mu_norm_tot; mu/max_fire]; % normalize by max fr
            end

            %% Plotting PSTHs for all units
            mu_avg = mean(mu_norm_tot,1);
            mu_std = std(mu_norm_tot,1);
            mu_sem = mu_std/sqrt(size(mu_avg,2));

            f2 = figure('Visible','off');
            x = 0:10:190;
            y = mu_avg;
            yu = mu_avg+mu_sem; 
            yl = mu_avg-mu_sem;
            fill([x fliplr(x)], [yu fliplr(yl)], [.9 .9 .9], 'linestyle', 'none');
            hold on;
            plot(x,y,'k');
            xline(0,'k--');
            dim = [.2 .5 .3 .3];
            str = 'N = 49';
            annotation('textbox',dim,'String',str,'FitBoxToText','on');
            xlabel('Time (10ms bins)');
            ylabel('Normalized Activity');
            xlim([-10 200]);
            ylim([0.4 0.9]);
            annotation('textarrow',[0.25 0.165],[0.2 0.2],'String','Motion Pulse Onset');
            title(sprintf('PSTH Across All Units: Epoch = %s, Pulse Speed = %d deg/s, Angle = %d',epoch,speed,angle));
            saveas(f2,sprintf('%s/figs/psth/allunits/p%2.2d-e%s-d%3.3d.png',stuff,speed,epoch,angle));
        end
    end
end

