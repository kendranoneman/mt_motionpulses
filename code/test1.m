clear;

a = [1; 1; 1; 1; 1; 1;2; 2; 2; 2; 2; 2;4; 4; 4; 4; 4; 4; 8; 8; 8; 8; 8; 8;16; 16; 16;16; 16; 16];
pulses = [1 2 4 8 16];

pulse_1 = find(a==1); 
pulse_2 = find(a==2);
pulse_4 = find(a==4);
pulse_8 = find(a==8);
pulse_16 = find(a==16);

for i=1:size(a,1)
    for j=1:size(pulses,2)
        pulse_{pulses(j)} = find(a==j);
    end
end