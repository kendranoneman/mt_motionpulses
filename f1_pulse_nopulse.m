clear all;

%% Load data
d = load('-mat','data/combinedMaestroSpkSortFEF(pa59pulsA).mat');

%% Finding names of good neurons

units_tot = [];
for y = 1:3372
   unit_vars1 = char(fieldnames(d.exp.dataMaestroPlx(y).units));
   unit_vals1 = struct2cell(d.exp.dataMaestroPlx(y).units); 
   
   unit_name_tot = [];
   for yy = 1:size(unit_vars1,1)
        unit_name = unit_vars1(yy,[5:8]);
        unit_name2 = {unit_name};
        unit_name_tot = [unit_name_tot; unit_name2];
   end
   units_tot = [units_tot; unit_name_tot];
end

[unique_units,~,idx] = unique(units_tot,'stable');
count_units = hist(idx,unique(idx));
max_num = max(count_units);

max_neurons = [];
for z = 1:size(count_units,2)
    if count_units(z) == max_num
        name = unique_units(z);
        max_neurons = [max_neurons; name];
    else
    end
end

%% Organization

%pulses = [1 2 4 8 16];
%angles = [0 90 180 270];
%epochs = ['f1'; 'f2'; 'p1'; 'p2'];

pulses = [16];
angles = [0 90 180 270];
epochs = ['f1'];

pulsef1 = 350; rangef1 = [pulsef1-150:pulsef1+250];
pulsef2 = 800;
pulsep1 = 1600;
pulsep2 = 2050; 

for g = 1:size(pulses,2)
    for gg = 1:size(angles,2)
        for ggg = 1:size(epochs,1)
            pulse = pulses(g);    % speed of motion pulse
            angle = angles(gg);
            epoch = epochs(ggg,:);

            values = struct2cell(d.exp.dataMaestroPlx);

            % file structures
            trName = char(values(1,:,:));
            trType = char(values(2,:,:));

            % empty matrices
            dir_pursuit_tot = [];
            speed_pursuit_tot = [];
            f1_tot = []; f2_tot = []; 
            p1_tot = []; p2_tot = [];
            pulse_duration_tot = []; pulse_speed_tot = [];
            pulse_speed_count = [];
            num_tot = [];

            % looping through each trial
            for i=1:size(trName,1)
                % trial name
                a1 = trName(i,:);
                num = str2num(a1(end-3:end)); num_tot = [num_tot; num];
                % trial type
                a2 = trType(i,:);
                dir_pursuit = str2num(a2(2:4)); dir_pursuit_tot = [dir_pursuit_tot; dir_pursuit];
                speed_pursuit = str2num(a2(8:9)); speed_pursuit_tot = [speed_pursuit_tot; speed_pursuit];
                f1 = str2num(a2(13:15));
                if isempty(f1)
                    f1 = 666;
                end
                f1_tot = [f1_tot; f1];
                f2 = str2num(a2(19:21)); 
                 if isempty(f2)
                    f2 = 666;
                end
                f2_tot = [f2_tot; f2];
                p1 = str2num(a2(25:27)); 
                if isempty(p1)
                    p1 = 666;
                end
                p1_tot = [p1_tot; p1];
                p2 = str2num(a2(31:33)); 
                if isempty(p2)
                    p2 = 666;
                end
                p2_tot = [p2_tot; p2];
                pulse_duration = str2num(a2(35:36)); pulse_duration_tot = [pulse_duration_tot; pulse_duration];
                pulse_speed = str2num(a2(37:38)); pulse_speed_tot = [pulse_speed_tot; pulse_speed];
            end

            speed_result = find(pulse_speed_tot==pulse);
            
            % f1-XXX-p1-XXX or f1-XXX-XXX-p2 vs. XXX-f2-XXX-p2
            if isequal(epoch,'f1') 
                y_pulse = find(f1_tot==angle & f2_tot==666 & p1_tot==angle & p2_tot==666);
                y2_pulse = find(f1_tot==angle & f2_tot==666 & p1_tot==666 & p2_tot==angle);
                
                n_pulse = find(f1_tot==666 & f2_tot==angle & p1_tot==666 & p2_tot==angle);
                
                p_result = union(y_pulse,y2_pulse,'sorted');
                n_result = n_pulse;
                
            % XXX-f2-XXX-p2 vs. f1-XXX-p1-XXX or f1-XXX-XXX-p2
            elseif isequal(epoch,'f2')
                y_pulse = find(f1_tot==666 & f2_tot==angle & p1_tot==666 & p2_tot==angle);
                
                n_pulse = find(f1_tot==angle & f2_tot==666 & p1_tot==angle & p2_tot==666);
                n2_pulse = find(f1_tot==angle & f2_tot==666 & p1_tot==666 & p2_tot==angle);
                
                p_result = y_pulse;
                n_result = union(n_pulse,n2_pulse,'sorted');
            
            % f1-XXX-p1-XXX vs. f1-XXX-XXX-p2 or XXX-f2-XXX-p2
            elseif isequal(epoch,'p1') 
                y_pulse = find(f1_tot==angle & f2_tot==666 & p1_tot==angle & p2_tot==666);
                
                n_pulse = find(f1_tot==angle & f2_tot==666 & p1_tot==666 & p2_tot==angle);
                n2_pulse = find(f1_tot==666 & f2_tot==angle & p1_tot==666 & p2_tot==angle);
                
                p_result = y_pulse;
                n_result = union(n_pulse,n2_pulse,'sorted');
                
            % f1-XXX-XXX-p2 or XXX-f2-XXX-p2 vs. f1-XXX-p1-XXX
            elseif isequal(epoch,'p2') 
                y_pulse = find(f1_tot==angle & f2_tot==666 & p1_tot==666 & p2_tot==angle);
                y2_pulse = find(f1_tot==666 & f2_tot==angle & p1_tot==666 & p2_tot==angle);
                
                n_pulse = find(f1_tot==angle & f2_tot==666 & p1_tot==angle & p2_tot==666);
               
                p_result = union(y_pulse,y2_pulse,'sorted');
                n_result = n_pulse;
                
             else
                fprintf('ew');
            end

            yes_pulse = intersect(speed_result,p_result);
            no_pulse = intersect(speed_result,n_result);


            %% Making D struct file for DimReduce
            units = values(10,:,:);

            D = struct([]);
            counter = 0;
            for ii = 1:size(yes_pulse,1)
                index = yes_pulse(ii);
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

                maxTime = 2851;
                N_tot = [];
                for l = 1:size(good_vals,1)
                    ca = good_vals(l); ca2 = [ca{:}]';
                    edges = 0 : 1 : (maxTime + 1);
                    [N,edges] = histcounts(ca2,edges);
                    N_tot(l,:) = N;
                end

                N_final = N_tot(:,:);
                D(ii).data = N_final;
                D(ii).condition = 'pulse';
                D(ii).epochStarts = [1, 350];
                %D(ii).epochColors = [0.87, 0.00,0.87];

            end

            start = size(yes_pulse,1);
            for jj = 1:size(no_pulse,1)
                spot = start+jj;
                index2 = no_pulse(jj);
                unit_vars = char(fieldnames(d.exp.dataMaestroPlx(index2).units));
                unit_vals = struct2cell(d.exp.dataMaestroPlx(index2).units);

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

                maxTime = 2851;
                N_tot = [];
                for l = 1:size(good_vals,1)
                    ca = good_vals(l); ca2 = [ca{:}]';
                    edges = 0 : 1 : (maxTime + 1);
                    [N,edges] = histcounts(ca2,edges);
                    N_tot(l,:) = N;
                end

                N_final = N_tot(:,:);
                D(spot).data = N_final;
                D(spot).condition = 'no pulse';
                D(spot).epochStarts = [1, 350];
                %D(spot).epochColors = [0.00,0.00,0.00];

            end

            % save struct D file
            save(sprintf('data/rawspiketrains-%s-p%d-d%d.mat',epoch,pulse,angle),'D');
        end
    end
end


