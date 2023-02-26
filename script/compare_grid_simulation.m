%% Load data

clear

time = datenum(2022,6,4,06,00,00);

str1 = datestr(time-1,"yyyymmddHHMM");
str2 = datestr(time,"yyyymmddHHMM");

mat = load("D:\MODEL_DATA_FSM\FSM_CURR\OUTPUT_GRID_0250\RESULTS_24h_opn\MODELDATA_" + str1 + "-" + str2 + "_FSM22.mat");
jul = load("D:\FSM_JULIA\" + str2 + "_output.mat");


%% Fix data

swe_mat = mat.swet.data;
swe_jul = jul.swe;

inan = swe_mat<0;

swe_mat(inan) = NaN;
swe_jul(inan) = NaN;


%% Plot results

fig = figure("Position",[294 151 1330 886]);

t = tiledlayout(2,2);

title(t,str2);

ax(1) = nexttile();
imagesc(flipud(swe_mat))
colorbar()
title("FSM Fortran/Matlab")
colormap('turbo');

ax(2) = nexttile();
imagesc(flipud(swe_jul))
colorbar()
title("FSM Julia")
colormap('turbo');

ax(3) = nexttile();
imagesc(flipud(swe_jul-swe_mat))
colorbar()
title("FSM Julia minus FSM Fortran/Matlab")
colormap('turbo');

nexttile()
plot(swe_mat(:),swe_jul(:),'.')
xlabel("Matlab/Fortran")
ylabel("Julia")

linkaxes(ax)

disp("Figures saved at: " + userpath)

saveas(fig,fullfile(userpath, time + "_swe.png"))

