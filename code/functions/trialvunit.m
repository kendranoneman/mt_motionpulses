function trialvunit(d,window,max_neurons,stuff)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
pulses = [1 2 4 8 16];
angles = [0 90 180 270];                                                                                         
epochs = ['f1'; 'f2'; 'p1'; 'p2'];

pulsef1 = 350; rangef1 = [pulsef1:(pulsef1+window)];
pulsef2 = 800; rangef2 = [pulsef2:(pulsef2+window)];
pulsep1 = 1600; rangep1 = [pulsep1:(pulsep1+window)];
pulsep2 = 2050; rangep2 = [pulsep2:(pulsep2+window)];

for p = 1:size(pulses,2)
    pulse = pulses(p);
    for a = 1:size(angles,2)
        angle = angles(a);
        for e = 1:size(epochs,1)
            epoch = epochs(e,:);
            fprintf('Pulse Speed = %d  Angle = %d  Epoch = %s\n',pulse,angle,epoch);
        
            values = struct2cell(d.exp.dataMaestroPlx);
            trName = char(values(1,:,:));
            trType = char(values(2,:,:));

            dir_pursuit_tot = []; speed_pursuit_tot = [];
            pulse_duration_tot = []; pulse_speed_tot = [];
            f1_tot = []; f2_tot = []; 
            p1_tot = []; p2_tot = [];

            % looping through each trial
            for i=1:size(trName,1)
                a2 = trType(i,:); % trial type
                dir_pursuit = str2num(a2(2:4)); dir_pursuit_tot = [dir_pursuit_tot; dir_pursuit];
                speed_pursuit = str2num(a2(8:9)); speed_pursuit_tot = [speed_pursuit_tot; speed_pursuit];

                f1 = str2num(a2(13:15)); % Fixation 1
                if isempty(f1)
                    f1 = 666;
                end
                f1_tot = [f1_tot; f1];
                f2 = str2num(a2(19:21)); % Fixation 2
                 if isempty(f2)
                    f2 = 666;
                end
                f2_tot = [f2_tot; f2];
                p1 = str2num(a2(25:27)); % Pursuit 1
                if isempty(p1)
                    p1 = 666;
                end
                p1_tot = [p1_tot; p1];
                p2 = str2num(a2(31:33)); % Pursuit 2
                if isempty(p2)
                    p2 = 666;
                end
                p2_tot = [p2_tot; p2];

                pulse_duration = str2num(a2(35:36)); pulse_duration_tot = [pulse_duration_tot; pulse_duration];
                pulse_speed = str2num(a2(37:38)); pulse_speed_tot = [pulse_speed_tot; pulse_speed];
            end

            % f1-XXX-p1-XXX or f1-XXX-XXX-p2 vs. XXX-f2-XXX-p2
            if isequal(epoch,'f1') 
                y_pulse = find(f1_tot==angle & f2_tot==666 & p1_tot==angle & p2_tot==666);
                y2_pulse = find(f1_tot==angle & f2_tot==666 & p1_tot==666 & p2_tot==angle);

                n_pulse = find(f1_tot==666 & f2_tot==angle & p1_tot==666 & p2_tot==angle);

                p_result = union(y_pulse,y2_pulse,'sorted');
                n_result = n_pulse;

                ind_range = rangef1; es_range = [1,150];

            % XXX-f2-XXX-p2 vs. f1-XXX-p1-XXX or f1-XXX-XXX-p2
            elseif isequal(epoch,'f2')
                y_pulse = find(f1_tot==666 & f2_tot==angle & p1_tot==666 & p2_tot==angle);

                n_pulse = find(f1_tot==angle & f2_tot==666 & p1_tot==angle & p2_tot==666);
                n2_pulse = find(f1_tot==angle & f2_tot==666 & p1_tot==666 & p2_tot==angle);

                p_result = y_pulse;
                n_result = union(n_pulse,n2_pulse,'sorted');

                ind_range = rangef2; es_range = [1,150];

            % f1-XXX-p1-XXX vs. f1-XXX-XXX-p2 or XXX-f2-XXX-p2
            elseif isequal(epoch,'p1') 
                y_pulse = find(f1_tot==angle & f2_tot==666 & p1_tot==angle & p2_tot==666);

                n_pulse = find(f1_tot==angle & f2_tot==666 & p1_tot==666 & p2_tot==angle);
                n2_pulse = find(f1_tot==666 & f2_tot==angle & p1_tot==666 & p2_tot==angle);

                p_result = y_pulse;
                n_result = union(n_pulse,n2_pulse,'sorted');

                ind_range = rangep1; es_range = [1,150];

            % f1-XXX-XXX-p2 or XXX-f2-XXX-p2 vs. f1-XXX-p1-XXX
            elseif isequal(epoch,'p2') 
                y_pulse = find(f1_tot==angle & f2_tot==666 & p1_tot==666 & p2_tot==angle);
                y2_pulse = find(f1_tot==666 & f2_tot==angle & p1_tot==666 & p2_tot==angle);

                n_pulse = find(f1_tot==angle & f2_tot==666 & p1_tot==angle & p2_tot==666);

                p_result = union(y_pulse,y2_pulse,'sorted');
                n_result = n_pulse;

                ind_range = rangep2; es_range = [1,150];
            end

            result_pulse = find(pulse_speed_tot==pulse);
            %result_angle = union(p_result,n_result,'sorted');
            result = intersect(result_pulse,p_result);
            
            
            D = struct([]); G = struct([]);
            N_counts = [];
            for ff = 1:size(result,1) %for each trial
                index = result(ff);
                unit_vars = char(fieldnames(d.exp.dataMaestroPlx(index).units));
                unit_vals = struct2cell(d.exp.dataMaestroPlx(index).units);

                good_vals = [];
                for q = 1:size(unit_vars,1)
                   n1 = unit_vars(q, [5:8]);
                   for qq = 1:size(max_neurons,1)
                      n2 = num2str(cell2mat(max_neurons(qq)));
                      if isequal(n1,n2)
                          val = unit_vals(q);
                          good_vals = [good_vals; val];
                      else
                      end
                   end
                end
                
                
                [N_tot,N_b20] = spiketrains_binning(good_vals,ind_range);
                N_final = N_tot(:,ind_range); % one trial, #units x 200
                
                G(ff).bin20 = N_b20;
                
                D(ff).data = N_final;
                D(ff).condition = sprintf('p%2.2d-e%s-d%3.3d',pulse,epoch,angle);
                D(ff).epochStarts = 1;
                
                N_sum = sum(N_final,2); %units x counts
                N_counts = [N_counts; N_sum'];
                
                if pulse==01
                    D(ff).epochColors = [76 0 153]/255;
                elseif pulse==02
                    D(ff).epochColors = [127 0 255]/255;
                elseif pulse==4
                    D(ff).epochColors = [153 51 255]/255;
                elseif pulse==8
                    D(ff).epochColors = [178 102 255]/255;
                elseif pulse==16
                    D(ff).epochColors = [204 153 255]/255;
                end
            end
            
            save(sprintf('%s/data/manifold/spikes/p%2.2d-e%s-d%3.3d.mat',stuff,pulse,epoch,angle),'G');
            save(sprintf('%s/data/manifold/counts/p%2.2d-e%s-d%3.3d.mat',stuff,pulse,epoch,angle),'N_counts');
            save(sprintf('%s/data/manifold/datahigh/p%2.2d-e%s-d%3.3d.mat',stuff,pulse,epoch,angle),'D'); 
        end
    end
end
end

