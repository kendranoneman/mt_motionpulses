clear all;
% Do NOT combine the no-pulse pursuit conditions
% different beast there since they already made an eye movement.

%% Load data
d = load('-mat','data/combinedMaestroSpkSortFEF(pa59pulsA).mat');
if not(isfolder('data/epochs/nopulse'))
    mkdir('data/epochs/nopulse')
end

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

pulses = [1 2 4 8 16];
%angles = [0 90 180 270];
epochs = ['f1'; 'f2'];

pulsef1 = 350; rangef1 = [pulsef1-150:pulsef1+250];
pulsef2 = 800; rangef2 = [pulsef2-150:pulsef2+250];
pulsep1 = 1600; rangep1 = [pulsep1-150:pulsep1+250];
pulsep2 = 2050; rangep2 = [pulsep2-150:pulsep2+250];

p = [0 4 8 12 16];
for g = 1:size(pulses,2)
    pp = p(g);
    for ggg = 1:size(epochs,1)
        pulse = pulses(g);    % speed of motion pulse
        epoch = epochs(ggg,:);

        num = pp+ggg; 
        update = ['Progress = ', num2str(num), '/20'];
        disp(update)
        
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
%             p1_tot = [p1_tot; p1];
%             p2 = str2num(a2(31:33)); 
%             if isempty(p2)
%                 p2 = 666;
%             end
%             p2_tot = [p2_tot; p2];
%             pulse_duration = str2num(a2(35:36)); pulse_duration_tot = [pulse_duration_tot; pulse_duration];
%             pulse_speed = str2num(a2(37:38)); pulse_speed_tot = [pulse_speed_tot; pulse_speed];
        end

        speed_result = find(pulse_speed_tot==pulse);

        % f1-XXX-p1-XXX or f1-XXX-XXX-p2 vs. XXX-f2-XXX-p2
        if isequal(epoch,'f1') 
            n_pulse_0 = find(f1_tot==666 & f2_tot==0 & p1_tot==666 & p2_tot==0);
            n_result_0 = n_pulse_0;
            
            n_pulse_90 = find(f1_tot==666 & f2_tot==90 & p1_tot==666 & p2_tot==90);
            n_result_90 = n_pulse_90;
            
            n_pulse_180 = find(f1_tot==666 & f2_tot==180 & p1_tot==666 & p2_tot==180);
            n_result_180 = n_pulse_180;
            
            n_pulse_270 = find(f1_tot==666 & f2_tot==270 & p1_tot==666 & p2_tot==270);
            n_result_270 = n_pulse_270;
            
            ind_range = rangef1;
            pulse_start = pulsef1;

        % XXX-f2-XXX-p2 vs. f1-XXX-p1-XXX or f1-XXX-XXX-p2
        elseif isequal(epoch,'f2')
            n_pulse_0 = find(f1_tot==0 & f2_tot==666 & p1_tot==0 & p2_tot==666);
            n2_pulse_0 = find(f1_tot==0 & f2_tot==666 & p1_tot==666 & p2_tot==0);
            n_result_0 = union(n_pulse_0,n2_pulse_0,'sorted');

            n_pulse_90 = find(f1_tot==90 & f2_tot==666 & p1_tot==90 & p2_tot==666);
            n2_pulse_90 = find(f1_tot==90 & f2_tot==666 & p1_tot==666 & p2_tot==90);
            n_result_90 = union(n_pulse_90,n2_pulse_90,'sorted');
            
            n_pulse_180 = find(f1_tot==180 & f2_tot==666 & p1_tot==180 & p2_tot==666);
            n2_pulse_180 = find(f1_tot==180 & f2_tot==666 & p1_tot==666 & p2_tot==180);
            n_result_180 = union(n_pulse_180,n2_pulse_180,'sorted');
                
            n_pulse_270 = find(f1_tot==270 & f2_tot==666 & p1_tot==270 & p2_tot==666);
            n2_pulse_270 = find(f1_tot==270 & f2_tot==666 & p1_tot==666 & p2_tot==270);
            n_result_270 = union(n_pulse_270,n2_pulse_270,'sorted');        
    
            ind_range = rangef2;
            pulse_start = pulsef2;

        % f1-XXX-p1-XXX vs. f1-XXX-XXX-p2 or XXX-f2-XXX-p2
%         elseif isequal(epoch,'p1') 
% 
%             n_pulse_0 = find(f1_tot==0 & f2_tot==666 & p1_tot==666 & p2_tot==0);
%             n2_pulse_0 = find(f1_tot==666 & f2_tot==0 & p1_tot==666 & p2_tot==0);
%             n_result_0 = union(n_pulse_0,n2_pulse_0,'sorted');
%             
%             n_pulse_90 = find(f1_tot==90 & f2_tot==666 & p1_tot==666 & p2_tot==90);
%             n2_pulse_90 = find(f1_tot==666 & f2_tot==90 & p1_tot==666 & p2_tot==90);
%             n_result_90 = union(n_pulse_90,n2_pulse_90,'sorted');
%             
%             n_pulse_180 = find(f1_tot==180 & f2_tot==666 & p1_tot==666 & p2_tot==180);
%             n2_pulse_180 = find(f1_tot==666 & f2_tot==180 & p1_tot==666 & p2_tot==180);
%             n_result_180 = union(n_pulse_180,n2_pulse_180,'sorted');
% 
%             n_pulse_270 = find(f1_tot==270 & f2_tot==666 & p1_tot==666 & p2_tot==270);
%             n2_pulse_270 = find(f1_tot==666 & f2_tot==270 & p1_tot==666 & p2_tot==270);
%             n_result_270 = union(n_pulse_270,n2_pulse_270,'sorted');
%             
%             ind_range = rangep1;
%             pulse_start = pulsep1;
% 
%         % f1-XXX-XXX-p2 or XXX-f2-XXX-p2 vs. f1-XXX-p1-XXX
%         elseif isequal(epoch,'p2') 
%             n_pulse_0 = find(f1_tot==0 & f2_tot==666 & p1_tot==0 & p2_tot==666);
%             n_result_0 = n_pulse_0;
%             
%             n_pulse_90 = find(f1_tot==90 & f2_tot==666 & p1_tot==90 & p2_tot==666);
%             n_result_90 = n_pulse_90;
%             
%             n_pulse_180 = find(f1_tot==180 & f2_tot==666 & p1_tot==180 & p2_tot==666);
%             n_result_180 = n_pulse_180;
%             
%             n_pulse_270 = find(f1_tot==270 & f2_tot==666 & p1_tot==270 & p2_tot==666);
%             n_result_270 = n_pulse_270;
%             
%             ind_range = rangep2;
%             pulse_start = pulsep2;

         else
            fprintf('ew');
        end

        no_pulse_0 = intersect(speed_result,n_result_0);
        no_pulse_90 = intersect(speed_result,n_result_90);
        no_pulse_180 = intersect(speed_result,n_result_180);
        no_pulse_270 = intersect(speed_result,n_result_270);


        %% Making D struct file for DimReduce
        units = values(10,:,:);

        D = struct([]);
        counter = 0;
        % degree = 0
        for ii = 1:size(no_pulse_0,1)
            index = no_pulse_0(ii);
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
                N_trim = N(ind_range);
                N_tot(l,:) = N_trim;
            end

            N_final = N_tot(:,:);
            D(ii).data = N_final;
            D(ii).condition = '000';
            D(ii).epochStarts = [1,150]; %pulse at 350
            D(ii).epochColors = [255, 102, 255; 153, 0, 153]/255; %magenta

        end

        % degree = 90
        start = size(no_pulse_0,1);
        for jj = 1:size(no_pulse_90,1)
            spot = start+jj;
            index = no_pulse_90(jj);
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
                N_trim = N(ind_range);
                N_tot(l,:) = N_trim;
            end

            N_final = N_tot(:,:);
            D(spot).data = N_final;
            D(spot).condition = '090';
            D(spot).epochStarts = [1,150]; %pulse at 350
            D(spot).epochColors = [102, 255, 102; 0, 153, 0]/255; %green
        
        end
        
        % degree = 180
        start2 = size(no_pulse_90,1);
        for kk = 1:size(no_pulse_180,1)
            spot2 = start+start2+kk;
            index = no_pulse_180(kk);
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
                N_trim = N(ind_range);
                N_tot(l,:) = N_trim;
            end

            N_final = N_tot(:,:);
            D(spot2).data = N_final;
            D(spot2).condition = '180';
            D(spot2).epochStarts = [1,150]; %pulse at 350
            D(spot2).epochColors = [255, 178, 102; 204, 102, 0]/255; %orange
        
        end
        
        % degree = 270
        start3 = size(no_pulse_180,1);
        for rr = 1:size(no_pulse_270,1)
            spot3 = start+start2+start3+rr;
            index = no_pulse_270(rr);
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
                N_trim = N(ind_range);
                N_tot(l,:) = N_trim;
            end

            N_final = N_tot(:,:);
            D(spot3).data = N_final;
            D(spot3).condition = '270';
            D(spot3).epochStarts = [1,150]; %pulse at 350
            D(spot3).epochColors = [102, 178, 255; 0, 76, 153]/255; %blue
        
        end
        
        % save struct D file
        save(sprintf('data/epochs/nopulse/rawspiketrains-%s-p%d.mat',epoch,pulse),'D');
    end
end


