%% Settings

clear

station = "MCH.BLS2";

station = replace(station,".","_");


%% compare julia against rerun using fortran code with txt input

fortran = readmatrix("..\fortran\output_64\output_" + station + "_run_from_julia.txt");
julia = readmatrix("..\fortran\output_julia\output_" + station + "_test.txt");

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

