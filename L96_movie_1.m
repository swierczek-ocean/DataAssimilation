function L96_movie_1(Array_SqEnKF,Array_4DVar,Array_EDA,Array_En4DVar,Array_True,Array_Obs,...
    color1,color2,color3,color4,color5,color6,dim1,dim2,dim3,coords,jump)

ACC_Colors

dim1obs = floor(dim1/2)+1;
dim2obs = floor(dim2/2)+1;
dim3obs = floor(dim3/2)+1;

Array1 = Array_SqEnKF([dim1,dim2,dim3],:);
Array2 = Array_4DVar([dim1,dim2,dim3],:);
Array3 = Array_EDA([dim1,dim2,dim3],:);
Array4 = Array_En4DVar([dim1,dim2,dim3],:);
Array5 = Array_True([dim1,dim2,dim3],:);
Array6 = Array_Obs([dim1obs,dim2obs,dim3obs],:);

exp_time = size(Array1,2);

dasz = 11;
obssz1 = 20;
obssz2 = 10;
truesz = 3.5;
lag = 2*jump;

figure()
set(gcf, 'Position', [25, 25, 1600, 900])
h1 = plot3(Array1(1,1),Array1(2,1),Array1(3,1),'.','Color',Color(:,color1),'MarkerSize',dasz);
hold on
h2 = plot3(Array2(1,1),Array2(2,1),Array2(3,1),'.','Color',Color(:,color2),'MarkerSize',dasz);
h3 = plot3(Array3(1,1),Array3(2,1),Array3(3,1),'.','Color',Color(:,color3),'MarkerSize',dasz);
h4 = plot3(Array4(1,1),Array4(2,1),Array4(3,1),'.','Color',Color(:,color4),'MarkerSize',dasz);
h5 = plot3(Array5(1,1),Array5(2,1),Array5(3,1),'.','Color',Color(:,color5),'MarkerSize',truesz);
h6 = plot3(Array6(1,1),Array6(2,1),Array6(3,1),'*','Color',Color(:,color6),'MarkerSize',obssz2);
plot3(Array6(1,1),Array6(2,1),Array6(3,1),'.','Color',Color(:,color6),'MarkerSize',obssz1);
xlabel('x')
ylabel('y')
zlabel('z')
axis(coords)
legend([h1(1),h2(1),h3(1),h4(1),h5(1),h6(1)],'SqEnKF','4DVar','EDA','En4DVar','true','obs')
hold off
set(gca, 'nextplot','replacechildren', 'Visible','on');

nFrames = 471;
vidObj = VideoWriter('L96_movie_1.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 15;
open(vidObj);
writeVideo(vidObj, getframe(gcf));

for ii=2:lag+1
    jj = floor(ii/jump)+1;
    h5 = plot3(Array5(1,1:ii),Array5(2,1:ii),Array5(3,1:ii),'Color',Color(:,color5),'LineWidth',truesz);
    hold on
    h2 = plot3(Array2(1,1:ii),Array2(2,1:ii),Array2(3,1:ii),'.','Color',Color(:,color2),'MarkerSize',dasz);
    h3 = plot3(Array3(1,1:ii),Array3(2,1:ii),Array3(3,1:ii),'.','Color',Color(:,color3),'MarkerSize',dasz);
    h4 = plot3(Array4(1,1:ii),Array4(2,1:ii),Array4(3,1:ii),'.','Color',Color(:,color4),'MarkerSize',dasz);
    h1 = plot3(Array1(1,1:ii),Array1(2,1:ii),Array1(3,1:ii),'.','Color',Color(:,color1),'MarkerSize',dasz);
    h6 = plot3(Array6(1,1:jj),Array6(2,1:jj),Array6(3,1:jj),'*','Color',Color(:,color6),'MarkerSize',obssz2);
    plot3(Array6(1,1:jj),Array6(2,1:jj),Array6(3,1:jj),'.','Color',Color(:,color6),'MarkerSize',obssz1);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis(coords)
    legend([h1(1),h2(1),h3(1),h4(1),h5(1),h6(1)],'SqEnKF','4DVar','EDA','En4DVar','true','obs')
    hold off
    drawnow()
    writeVideo(vidObj, getframe(gcf));
end

for ii=lag+2:exp_time
    jj = floor(ii/jump)+1;
    kk = floor((ii-lag)/jump)+1;
    ll = ii-lag;
    h5 = plot3(Array5(1,ll:ii),Array5(2,ll:ii),Array5(3,ll:ii),'Color',Color(:,color5),'LineWidth',truesz);
    hold on
    h2 = plot3(Array2(1,ll:ii),Array2(2,ll:ii),Array2(3,ll:ii),'.','Color',Color(:,color2),'MarkerSize',dasz);
    h3 = plot3(Array3(1,ll:ii),Array3(2,ll:ii),Array3(3,ll:ii),'.','Color',Color(:,color3),'MarkerSize',dasz);
    h4 = plot3(Array4(1,ll:ii),Array4(2,ll:ii),Array4(3,ll:ii),'.','Color',Color(:,color4),'MarkerSize',dasz);
    h1 = plot3(Array1(1,ll:ii),Array1(2,ll:ii),Array1(3,ll:ii),'.','Color',Color(:,color1),'MarkerSize',dasz);
    h6 = plot3(Array6(1,kk:jj),Array6(2,kk:jj),Array6(3,kk:jj),'*','Color',Color(:,color6),'MarkerSize',obssz2);
    plot3(Array6(1,kk:jj),Array6(2,kk:jj),Array6(3,kk:jj),'.','Color',Color(:,color6),'MarkerSize',obssz1);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis(coords)
    legend([h1(1),h2(1),h3(1),h4(1),h5(1),h6(1)],'SqEnKF','4DVar','EDA','En4DVar','true','obs')
    hold off
    drawnow()
    writeVideo(vidObj, getframe(gcf));
end

close(vidObj);
end

