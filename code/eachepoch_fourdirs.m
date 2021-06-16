clear vars;

%% Load data
stuff = '/Users/kendranoneman/Projects/mayo';
d = load('-mat',sprintf('%s/data/combinedMaestroSpkSortFEF(pa59pulsA).mat',stuff));
if not(isfolder(sprintf('%s/data/eachepoch/fourdirs',stuff)))
    mkdir(sprintf('%s/data/eachepoch/fourdirs',stuff));
end

%% Finding names of good neurons
max_neurons = good_neurons(d);

%% Organization

pulses = [1 2 4 8 16];
%angles = [0 90 180 270];
epochs = ['f1'; 'f2'; 'p1'; 'p2'];

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
        dir_pursuit_tot = []; speed_pursuit_tot = [];
        f1_tot = []; f2_tot = []; p1_tot = []; p2_tot = [];
        pulse_duration_tot = []; pulse_speed_tot = [];
        pulse_speed_count = []; num_tot = [];

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
            % yes motion pulse
            y_pulse_0 = find(f1_tot==0 & f2_tot==666 & p1_tot==0 & p2_tot==666);
            y2_pulse_0 = find(f1_tot==0 & f2_tot==666 & p1_tot==666 & p2_tot==0);
            p_result_0 = union(y_pulse_0,y2_pulse_0,'sorted');
            
            y_pulse_90 = find(f1_tot==90 & f2_tot==666 & p1_tot==90 & p2_tot==666);
            y2_pulse_90 = find(f1_tot==90 & f2_tot==666 & p1_tot==666 & p2_tot==90);
            p_result_90 = union(y_pulse_90,y2_pulse_90,'sorted');
            
            y_pulse_180 = find(f1_tot==180 & f2_tot==666 & p1_tot==180 & p2_tot==666);
            y2_pulse_180 = find(f1_tot==180 & f2_tot==666 & p1_tot==666 & p2_tot==180);
            p_result_180 = union(y_pulse_180,y2_pulse_180,'sorted');
            
            y_pulse_270 = find(f1_tot==270 & f2_tot==666 & p1_tot==270 & p2_tot==666);
            y2_pulse_270 = find(f1_tot==270 & f2_tot==666 & p1_tot==666 & p2_tot==270);
            p_result_270 = union(y_pulse_270,y2_pulse_270,'sorted');
            
            % no motion pulse
            n_pulse_0 = find(f1_tot==666 & f2_tot==0 & p1_tot==666 & p2_tot==0);
            n_result_0 = n_pulse_0;
            
            n_pulse_90 = find(f1_tot==666 & f2_tot==90 & p1_tot==666 & p2_tot==90);
            n_result_90 = n_pulse_90;
            
            n_pulse_180 = find(f1_tot==666 & f2_tot==180 & p1_tot==666 & p2_tot==180);
            n_result_180 = n_pulse_180;
            
            n_pulse_270 = find(f1_tot==666 & f2_tot==270 & p1_tot==666 & p2_tot==270);
            n_result_270 = n_pulse_270;
            
            ind_range = rangef1;

        % XXX-f2-XXX-p2 vs. f1-XXX-p1-XXX or f1-XXX-XXX-p2
        elseif isequal(epoch,'f2')
            % yes motion pulse
            y_pulse_0 = find(f1_tot==666 & f2_tot==0 & p1_tot==666 & p2_tot==0);
            p_result_0 = y_pulse_0;
            
            y_pulse_90 = find(f1_tot==666 & f2_tot==90 & p1_tot==666 & p2_tot==90);
            p_result_90 = y_pulse_90;
            
            y_pulse_180 = find(f1_tot==666 & f2_tot==180 & p1_tot==666 & p2_tot==180);
            p_result_180 = y_pulse_180;
            
            y_pulse_270 = find(f1_tot==666 & f2_tot==270 & p1_tot==666 & p2_tot==270);
            p_result_270 = y_pulse_270;
            
            % no motion pulse
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

        % f1-XXX-p1-XXX vs. f1-XXX-XXX-p2 or XXX-f2-XXX-p2
        elseif isequal(epoch,'p1') 
            y_pulse_0 = find(f1_tot==0 & f2_tot==666 & p1_tot==0 & p2_tot==666);
            p_result_0 = y_pulse_0;
            
            y_pulse_90 = find(f1_tot==90 & f2_tot==666 & p1_tot==90 & p2_tot==666);
            p_result_90 = y_pulse_90;
            
            y_pulse_180 = find(f1_tot==180 & f2_tot==666 & p1_tot==180 & p2_tot==666);
            p_result_180 = y_pulse_180;
            
            y_pulse_270 = find(f1_tot==270 & f2_tot==666 & p1_tot==270 & p2_tot==666);
            p_result_270 = y_pulse_270;

            ind_range = rangep1;

        % f1-XXX-XXX-p2 or XXX-f2-XXX-p2 vs. f1-XXX-p1-XXX
        elseif isequal(epoch,'p2') 
            y_pulse_0 = find(f1_tot==0 & f2_tot==666 & p1_tot==666 & p2_tot==0);
            y2_pulse_0 = find(f1_tot==666 & f2_tot==0 & p1_tot==666 & p2_tot==0);
            p_result_0 = union(y_pulse_0,y2_pulse_0,'sorted');

            y_pulse_90 = find(f1_tot==90 & f2_tot==666 & p1_tot==666 & p2_tot==90);
            y2_pulse_90 = find(f1_tot==666 & f2_tot==90 & p1_tot==666 & p2_tot==90);
            p_result_90 = union(y_pulse_90,y2_pulse_90,'sorted');
            
            y_pulse_180 = find(f1_tot==180 & f2_tot==666 & p1_tot==666 & p2_tot==180);
            y2_pulse_180 = find(f1_tot==666 & f2_tot==180 & p1_tot==666 & p2_tot==180);
            p_result_180 = union(y_pulse_180,y2_pulse_180,'sorted');
            
            y_pulse_270 = find(f1_tot==270 & f2_tot==666 & p1_tot==666 & p2_tot==270);
            y2_pulse_270 = find(f1_tot==666 & f2_tot==270 & p1_tot==666 & p2_tot==270);
            p_result_270 = union(y_pulse_270,y2_pulse_270,'sorted');
            
            ind_range = rangep2;
        end

        yes_pulse_0 = intersect(speed_result,p_result_0);
        yes_pulse_90 = intersect(speed_result,p_result_90);
        yes_pulse_180 = intersect(speed_result,p_result_180);
        yes_pulse_270 = intersect(speed_result,p_result_270);
        
        no_pulse_0 = intersect(speed_result,n_result_0);
        no_pulse_90 = intersect(speed_result,n_result_90);
        no_pulse_180 = intersect(speed_result,n_result_180);
        no_pulse_270 = intersect(speed_result,n_result_270);
        
        no_pulse_all = union(union(union(no_pulse_0,no_pulse_90,'sorted'),no_pulse_180,'sorted'),no_pulse_270,'sorted');

        %% Making D struct file for DimReduce
        units = values(10,:,:);

        D = struct([]);
        start = 0;
        % degree = 0
        for jj = 1:size(yes_pulse_0,1)
            spot = start+jj;
            index = yes_pulse_0(jj);
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
            
            N_tot = spiketrains_binning(good_vals);

            N_final = N_tot(:,ind_range);
            D(spot).data = N_final;
            D(spot).condition = '000';
            D(spot).epochStarts = [1,150]; %pulse at 350
            D(spot).epochColors = [153, 0, 153; 255, 102, 255]/255; %magenta
        end
        
        % degree = 90
        for jj = 1:size(yes_pulse_90,1)
            spot2 = spot+jj;
            index = yes_pulse_90(jj);
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

            N_tot = spiketrains_binning(good_vals);

            N_final = N_tot(:,ind_range);
            D(spot2).data = N_final;
            D(spot2).condition = '090';
            D(spot2).epochStarts = [1,150]; %pulse at 350
            D(spot2).epochColors = [0, 153, 0; 102, 255, 102]/255; %green
        end
        
        % degree = 180
        for jj = 1:size(yes_pulse_180,1)
            spot3 = spot2+jj;
            index = yes_pulse_180(jj);
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

            N_tot = spiketrains_binning(good_vals);

            N_final = N_tot(:,ind_range);
            D(spot3).data = N_final;
            D(spot3).condition = '180';
            D(spot3).epochStarts = [1,150]; %pulse at 350
            D(spot3).epochColors = [204, 102, 0; 255, 178, 102]/255; %orange
        end
        
        % degree = 270
        for jj = 1:size(yes_pulse_270,1)
            spot4 = spot3+jj;
            index = yes_pulse_270(jj);
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

            N_tot = spiketrains_binning(good_vals);

            N_final = N_tot(:,ind_range);
            D(spot4).data = N_final;
            D(spot4).condition = '270';
            D(spot4).epochStarts = [1,150]; 
            D(spot4).epochColors = [0, 76, 153; 102, 178, 255]/255; %blue
        end
        
        if isequal(epoch,'f1') || isequal(epoch, 'f2')
            % no pulse
            for jj = 1:size(no_pulse_all,1)
                spot5 = spot4+jj;
                index = no_pulse_all(jj);
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

                N_tot = spiketrains_binning(good_vals);

                N_final = N_tot(:,ind_range);
                D(spot5).data = N_final;
                D(spot5).condition = 'no pulse';
                D(spot5).epochStarts = [1,150]; 
                D(spot5).epochColors = [0, 0, 0; 0, 0, 0]/255; %black
            end
        end
        
        % save struct D file
        save(sprintf('%s/data/eachepoch/fourdirs/rawspiketrains-%s-p%d.mat',stuff,epoch,pulse),'D');
    end
end


