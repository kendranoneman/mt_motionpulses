function [max_neurons] = good_neurons(d)
%% Finding names of good neurons
units_tot = [];
for y = 1:3372 %number of trials
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
end