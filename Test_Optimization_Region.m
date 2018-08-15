colormap(jet)
r = [2:0.2:8];             % localization radius
alpha = [0.0:0.01:0.16];   % ensemble inflation parameter
surf(alpha,r,sqrt(error_list_SqEnKF))
axis([0 0.16 2 8 0.3 1])
caxis([0.6 1])
