%% compare against jim operational runs (see branch jan-translate-input in jim_operational)

new = readmatrix("..\fortran\data\output_5wj_32.txt");
old = readmatrix("..\fortran\data\results_jim_operational_5wj.txt");

new_time = datenum(new(:,1),new(:,2),new(:,3),new(:,4),0,0);
old_time = datenum(old(:,1),old(:,2),old(:,3),old(:,4),0,0);

[~,inew,iold] = intersect(new_time,old_time);

figure
plot(old_time(iold),old(iold,5),'.r')
hold on
plot(new_time(inew),new(inew,5),'b')
title("Ds")

disp("Max diff in Ds :" + max(abs(new(inew,5)-old(iold,5))))

figure
plot(old_time(iold),old(iold,6),'.r')
hold on
plot(new_time(inew),new(inew,7),'b')
title("SWE")

disp("Max diff in SWE :" + max(abs(old(iold,6)-new(inew,7))))

figure
plot(old_time(iold),old(iold,7),'.r')
hold on
plot(new_time(inew),new(inew,8),'b')
title("Tsrf")

disp("Max diff in Tsrf :" + max(abs(old(iold,7)-new(inew,8))))

clear


%% verify restart of fsm...

complete = readmatrix("..\fortran\data\output_5wj_32.txt");
first_half = readmatrix("..\fortran\data\output_first_5wj_32.txt");
second_half = readmatrix("..\fortran\data\output_last_5wj_32.txt");
merged = [first_half; second_half];

disp("Check that all are the same: " + all(complete(:) == merged(:)))

clear

%% compare 32 and 64 bit versions

complete_32 = readmatrix("..\fortran\data\output_5wj_32.txt");
complete_64 = readmatrix("..\fortran\data\output_5wj_64.txt");

disp("Check that all are the same: " + all(complete_32(:) == complete_64(:)))

diff_hs    = complete_32(:,5) - complete_64(:,5);
diff_fsnow = complete_32(:,6) - complete_64(:,6);
diff_swe   = complete_32(:,7) - complete_64(:,7);
diff_tsrf  = complete_32(:,8) - complete_64(:,8);

figure
plot(diff_hs)
ylabel("diff hs")

figure
plot(diff_fsnow)
ylabel("diff fsnow")

figure
plot(diff_swe)
ylabel("diff swe")

figure
plot(diff_tsrf)
ylabel("diff tsrf")




