colormap(jet)
r2 = [0.5:0.25:8];      
alpha2 = [0:0.01:0.25];  
surf(alpha2,r2,sqrt(error_list_En4DVar))
axis([0 0.25 0.5 8 0 5])
caxis([0 5])
