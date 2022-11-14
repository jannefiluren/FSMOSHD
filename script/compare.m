%% compare against jim operational runs (see branch jan-translate-input in jim_operational)

new = readmatrix("..\fortran\data\output_5wj.txt");
old = readmatrix("..\fortran\data\results_jim_operational_5wj.txt");

new_time = datenum(new(:,1),new(:,2),new(:,3),new(:,4),0,0);
old_time = datenum(old(:,1),old(:,2),old(:,3),old(:,4),0,0);

figure
plot(old_time,old(:,5),'.r')
hold on
plot(new_time,new(:,5))
title("Ds")

figure
plot(new_time,new(:,7))
hold on
plot(old_time,old(:,6),'*r')
title("SWE")

figure
plot(new_time,new(:,8))
hold on
plot(old_time,old(:,7),'*r')
title("Tsrf")

clear

%% verify restart of fsm...

complete = readmatrix("..\fortran\data\output_5wj.txt");
first_half = readmatrix("..\fortran\data\output_first_5wj.txt");
second_half = readmatrix("..\fortran\data\output_last_5wj.txt");
merged = [first_half; second_half];

disp("Check that all are the same: " +all(complete(:) == merged(:))





