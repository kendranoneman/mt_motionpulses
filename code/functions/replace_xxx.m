function [f1_tot,f2_tot,p1_tot,p2_tot] = replace_xxx(a2,f1_tot,f2_tot,p1_tot,p2_tot)
%% Replacing XXX (no motion pulse) with 666 to keep arrays full of doubles only
% Chose 666, since I didn't want it to include 0,1,2,7,8, or 9 and wanted
% it to be the same length as the others (000, 090, 180, 270).
% Input = a2 --> trial type column from data file
% Output = f1_tot, f2_tot, p1_tot, p2_tot --> finished arrays with 666

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
end
