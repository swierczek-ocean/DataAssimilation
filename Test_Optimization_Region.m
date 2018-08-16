colormap(jet)
r = [1:0.2:8];            
alpha = [0.0:0.01:0.35]; 
surf(alpha,r,error_list_SqEnKF)
axis([0 0.35 1 8 0.3 1])
caxis([0.3 1])
