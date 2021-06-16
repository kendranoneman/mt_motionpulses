clear vars;

%% Load data
stuff = '/Users/kendranoneman/Projects/mayo';
d = load('-mat',sprintf('%s/data/combinedMaestroSpkSortFEF(pa59pulsA).mat',stuff));
if not(isfolder(sprintf('%s/data/eachepoch/pulsespeeds',stuff)))
    mkdir(sprintf('%s/data/eachepoch/pulsespeeds',stuff));
end

%% Finding names of good neurons
max_neurons = good_neurons(d);

%% Organization

pulses = [1 2 4 8 16];
angles = [0 90 180 270];
epochs = ['f1'; 'f2'; 'p1'; 'p2'];

pulsef1 = 350; rangef1 = [pulsef1-150:pulsef1+250];
pulsef2 = 800; rangef2 = [pulsef2-150:pulsef2+250];
pulsep1 = 1600; rangep1 = [pulsep1-150:pulsep1+250];
pulsep2 = 2050; rangep2 = [pulsep2-150:pulsep2+250];

%% Looping through angles and epochs
p = [0 4 8 12];
for gg = 1:size(angles,2)
    pp = p(gg);
    for ggg = 1:size(epochs,1)
        % speed of motion pulse
        angle = angles(gg);
        epoch = epochs(ggg,:);

        num = pp+ggg; 
        update = ['Progress = ', num2str(num), '/16'];
        disp(update)

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

        yes_pulse = p_result; no_pulse = n_result;
        
        pulse_1 = find(pulse_speed_tot==1); 
        pulse_2 = find(pulse_speed_tot==2);
        pulse_4 = find(pulse_speed_tot==4);
        pulse_8 = find(pulse_speed_tot==8);
        pulse_16 = find(pulse_speed_tot==16);

        yes_1 = intersect(yes_pulse,pulse_1); no_1 = intersect(no_pulse,pulse_1);
        yes_2 = intersect(yes_pulse,pulse_2); no_2 = intersect(no_pulse,pulse_2);
        yes_4 = intersect(yes_pulse,pulse_4); no_4 = intersect(no_pulse,pulse_4);
        yes_8 = intersect(yes_pulse,pulse_8); no_8 = intersect(no_pulse,pulse_8);
        yes_16 = intersect(yes_pulse,pulse_16); no_16 = intersect(no_pulse, pulse_16);

        no_pulse_tot = union(union(union(union(no_1,no_2,'sorted'),no_4,'sorted'),no_8,'sorted'),no_16,'sorted');

        %% Making D struct file for DimReduce
        units = values(10,:,:);

        D = struct([]);

        % Pulse occurs
        % speed = 1;
        start = 0;
        for aa = 1:size(yes_1,1)
            spot1 = start+aa;
            index1 = yes_1(aa);
            unit_vars = char(fieldnames(d.exp.dataMaestroPlx(index1).units));
            unit_vals = struct2cell(d.exp.dataMaestroPlx(index1).units);

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
            D(spot1).data = N_final;
            D(spot1).condition = '1';
            D(spot1).epochStarts = es_range;
            D(spot1).epochColors = [153, 0, 0; 76, 0, 153]/255; 
        end

        % speed = 2;
        for bb = 1:size(yes_2,1)
            spot2 = spot1+bb;
            index2 = yes_2(bb);
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

            N_tot = spiketrains_binning(good_vals);

            N_final = N_tot(:,ind_range);
            D(spot2).data = N_final;
            D(spot2).condition = '2';
            D(spot2).epochStarts = es_range;
            D(spot2).epochColors = [255, 0, 0; 127, 0, 255]/255; 
        end

        % speed = 4;
        for cc = 1:size(yes_4,1)
            spot4 = spot2+cc;
            index4 = yes_4(cc);
            unit_vars = char(fieldnames(d.exp.dataMaestroPlx(index4).units));
            unit_vals = struct2cell(d.exp.dataMaestroPlx(index4).units);

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
            D(spot4).condition = '4';
            D(spot4).epochStarts = es_range;
            D(spot4).epochColors = [233, 51, 51; 153, 51, 255]/255; 
        end

        % speed = 8;
        for dd = 1:size(yes_8,1)
            spot8 = spot4+dd;
            index8 = yes_8(dd);
            unit_vars = char(fieldnames(d.exp.dataMaestroPlx(index8).units));
            unit_vals = struct2cell(d.exp.dataMaestroPlx(index8).units);

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
            D(spot8).data = N_final;
            D(spot8).condition = '8';
            D(spot8).epochStarts = es_range;
            D(spot8).epochColors = [233, 51, 51; 153, 51, 255]/255; 
        end

        % speed = 16;
        for ee = 1:size(yes_16,1)
            spot16 = spot8+ee;
            index16 = yes_16(ee);
            unit_vars = char(fieldnames(d.exp.dataMaestroPlx(index16).units));
            unit_vals = struct2cell(d.exp.dataMaestroPlx(index16).units);

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
            D(spot16).data = N_final;
            D(spot16).condition = '16';
            D(spot16).epochStarts = es_range;
            D(spot16).epochColors = [255, 153, 153; 204, 153, 255]/255; 

        end

        % no pulse
        for ff = 1:size(no_pulse_tot,1)
            spotno = spot16+ff;
            indexno = no_pulse_tot(ff);
            unit_vars = char(fieldnames(d.exp.dataMaestroPlx(indexno).units));
            unit_vals = struct2cell(d.exp.dataMaestroPlx(indexno).units);

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
            D(spotno).data = N_final;
            D(spotno).condition = 'no pulse';
            D(spotno).epochStarts = es_range;
            D(spotno).epochColors = [0, 0, 0; 0, 0, 0]/255; %black
        end

        % save struct D file
        save(sprintf('%s/data/eachepoch/pulsespeeds/rawspiketrains-allp-%s-d%d.mat',stuff,epoch,angle),'D');
    end
end


