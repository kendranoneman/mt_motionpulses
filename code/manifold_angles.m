clear all;
folder = '/Users/kendranoneman/Projects/mayo/figs/manifold/stackedspeeds';

epochs = ['f1'; 'f2'; 'p1'; 'p2']; 
for epo = 1:size(epochs,1)  
    A = struct([]);
    for i=1:2
        epoch = epochs(epo,:); 
        load('-mat', sprintf('%s/xyz/xyz-e%s.mat',folder,epoch));

        format long;
        if i == 1
            v1 = double(XYZ(1).v1)'; v1 = v1/norm(v1);
            v2 = double(XYZ(2).v1)'; v2 = v2/norm(v2);
            v4 = double(XYZ(3).v1)'; v4 = v4/norm(v4);
            v8 = double(XYZ(4).v1)'; v8 = v8/norm(v8);
            v16 = double(XYZ(5).v1)'; v16 = v16/norm(v16);
        end
        if i == 2
            v1 = double(XYZ(1).v2)'; v1 = v1/norm(v1);
            v2 = double(XYZ(2).v2)'; v2 = v2/norm(v2);
            v4 = double(XYZ(3).v2)'; v4 = v4/norm(v4);
            v8 = double(XYZ(4).v2)'; v8 = v8/norm(v8);
            v16 = double(XYZ(5).v2)'; v16 = v16/norm(v16);
        end

        A(i).ang1_2 = atan2(norm(cross(v1,v2)),dot(v1,v2))*(180/pi);
        A(i).ang1_4 = atan2(norm(cross(v1,v4)),dot(v1,v4))*(180/pi);
        A(i).ang1_8 = atan2(norm(cross(v1,v8)),dot(v1,v8))*(180/pi);
        A(i).ang1_16 = atan2(norm(cross(v1,v16)),dot(v1,v16))*(180/pi);
        A(i).ang2_4 = atan2(norm(cross(v2,v4)),dot(v2,v4))*(180/pi);
        A(i).ang2_8 = atan2(norm(cross(v2,v8)),dot(v2,v8))*(180/pi);
        A(i).ang2_16 = atan2(norm(cross(v2,v16)),dot(v2,v16))*(180/pi);
        A(i).ang4_8 = atan2(norm(cross(v4,v8)),dot(v4,v8))*(180/pi);
        A(i).ang4_16 = atan2(norm(cross(v4,v16)),dot(v4,v16))*(180/pi);
        A(i).ang8_16 = atan2(norm(cross(v8,v16)),dot(v8,v16))*(180/pi);
        
        save(sprintf('%s/angles/e%s.mat',folder,epoch), 'A');
    end
end

%% Histogram of angles across epochs
f1 = figure;
ef1 = load('-mat',sprintf('%s/angles/ef1.mat',folder));
ef2 = load('-mat',sprintf('%s/angles/ef2.mat',folder));
ep1 = load('-mat',sprintf('%s/angles/ep1.mat',folder));
ep2 = load('-mat',sprintf('%s/angles/ep2.mat',folder));
angles_f1 = struct2table(ef1.A); a_f1 = angles_f1{1,1:10}';
angles_f2 = struct2table(ef2.A); a_f2 = angles_f2{1,1:10}';
angles_p1 = struct2table(ep1.A); a_p1 = angles_p1{1,1:10}';
angles_p2 = struct2table(ep2.A); a_p2 = angles_p2{1,1:10}';
aa = [a_f1 a_f2 a_p1 a_p2];
