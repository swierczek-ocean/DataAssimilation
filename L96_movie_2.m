function L96_movie_2(Array_True,color1,dim1,dim2,dim3,coords)

ACC_Colors

Array = Array_True([dim1,dim2,dim3],:);

exp_time = size(Array,2);

truesz = 3.5;

figure()
set(gcf, 'Position', [25, 25, 1600, 900])
plot3(Array(1,1),Array(2,1),Array(3,1),'.','Color',Color(:,color1),'MarkerSize',truesz);
xlabel('x')
ylabel('y')
zlabel('z')
axis(coords)
set(gca, 'nextplot','replacechildren', 'Visible','on');

nFrames = 471;
vidObj = VideoWriter('L96_movie_2.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 25;
open(vidObj);
writeVideo(vidObj, getframe(gcf));

for ii=2:exp_time
    plot3(Array(1,1:ii),Array(2,1:ii),Array(3,1:ii),'Color',Color(:,color1),'LineWidth',truesz);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis(coords)
    drawnow()
    writeVideo(vidObj, getframe(gcf));
end

close(vidObj);
end

