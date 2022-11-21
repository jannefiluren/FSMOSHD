%% compare against jim operational runs (see branch jan-translate-input in jim_operational)

new = readmatrix("..\fortran\output_32\output_SLF_5WJ.txt");
old = readmatrix("..\fortran\output_oshd\results_SLF_5WJ.txt");

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

complete = readmatrix("..\fortran\output_32\output_SLF_5WJ.txt");
first_half = readmatrix("..\fortran\output_32\output_first_SLF_5WJ.txt");
second_half = readmatrix("..\fortran\output_32\output_last_SLF_5WJ.txt");
merged = [first_half; second_half];

disp("Check that all are the same: " + all(complete(:) == merged(:)))

clear

%% compare 32 and 64 bit versions

complete_32 = readmatrix("..\fortran\output_32\output_SLF_5WJ.txt");
complete_64 = readmatrix("..\fortran\output_64\output_SLF_5WJ.txt");

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


%% compare julia against rerun using fortran code with txt input

fortran = readmatrix("..\fortran\output_64\output_SLF_5WJ.txt");
julia = readmatrix("..\fortran\output_julia\output_SLF_5WJ.txt");

fortran_time = datenum(fortran(:,1),fortran(:,2),fortran(:,3),fortran(:,4),0,0);
julia_time = datenum(julia(:,1),julia(:,2),julia(:,3),julia(:,4),0,0);

[~,ifortran,ijulia] = intersect(fortran_time,julia_time);

figure
plot(julia_time(ijulia),julia(ijulia,5),'.r')
hold on
plot(fortran_time(ifortran),fortran(ifortran,5),'b')
legend("julia","fortran")
title("Ds")

disp("Max diff in Ds :" + max(abs(julia(ijulia,5)-fortran(ifortran,5))))

figure
plot(julia_time(ijulia),julia(ijulia,7),'.r')
hold on
plot(fortran_time(ifortran),fortran(ifortran,7),'b')
legend("julia","fortran")
title("SWE")

disp("Max diff in SWE :" + max(abs(julia(ijulia,7)-fortran(ifortran,7))))

figure
plot(julia_time(ijulia),julia(ijulia,8),'.r')
hold on
plot(fortran_time(ifortran),fortran(ifortran,8),'b')
legend("julia","fortran")
title("Tsrf")

disp("Max diff in Tsrf :" + max(abs(julia(ijulia,8)-fortran(ifortran,8))))

figure
plot(julia_time(ijulia),julia(ijulia,9),'.r')
hold on
plot(fortran_time(ifortran),fortran(ifortran,9),'b')
legend("julia","fortran")
title("Nsnow")

disp("Max diff in Nsnow :" + max(abs(julia(ijulia,9)-fortran(ifortran,9))))

