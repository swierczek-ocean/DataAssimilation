%colormap(jet)
figure()
r = [1:0.2:8];            
alpha = [0.0:0.01:0.35]; 
surf(alpha,r,error_list_SqEnKF)
axis([0 0.35 1 8 0.3 1])
title('SqEnKF')
caxis([0.3 1])
% 
% 
% figure()
% r = [1:0.2:8];            
% alpha = [0.0:0.01:0.35]; 
% surf(alpha,r,error_list_EDA)
% axis([0 0.35 1 8 0.3 1.5])
% title('EDA')
% caxis([0.3 1.5])


figure()
r = [1:0.2:8];            
alpha = [0.0:0.01:0.35]; 
surf(alpha,r,error_list_En4DVar)
axis([0 0.35 1 8 0.3 6])
title('En4DVar')
caxis([0.3 6])


figure()
r = [1:0.2:8];            
alpha = [0.0:0.01:0.35]; 
surf(alpha,r,error_list_4DVar)
axis([0 0.35 1 8 0.3 2])
title('4DVar')
caxis([0.3 2])




