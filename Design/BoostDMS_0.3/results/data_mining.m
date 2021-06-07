    clear;
close all;
load BoostDMS_lastiteration_problem1.mat;

protoype = [0.05257; 0.04814; 0.1624; 0.0832; 0.1175; 0.027;0.3491; pi/2];
[ d, ix ] = min( vecnorm( Plist-protoype, 1 ) );

indexes = Flist(3,:) < 30;

to_study_points = Plist(:,not(indexes));
to_study_points_f_values = Flist(:,not(indexes));

for i = 1:length(Plist)
    i;
    Plist(:,i)
    SixRSS (Plist(:,i));
    input('')
    close all
end



Plist(ix)

NewFlist = Flist(:,indexes);
NewPlist = Plist(:,indexes);
figure
hold on
scatter3(NewFlist(1,:),NewFlist(2,:),NewFlist(3,:))
xlabel('Vol')
ylabel('GTSI')
zlabel('GRSI')