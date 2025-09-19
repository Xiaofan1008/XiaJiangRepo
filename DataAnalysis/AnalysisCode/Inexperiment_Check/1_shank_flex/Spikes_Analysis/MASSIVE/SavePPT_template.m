%%overwites ppt with name in presntationPath!!!
import mlreportgen.ppt.* %need this to import ppt save format
TemplateFile = 'C:\Users\smei0006\Documents\myRasterTemplate.pptx'; %template where you can alter slide master and selection pane names layout etc.
presentationPath = 'testRaster.pptx'; %saving file
presentationObj = Presentation(presentationPath,TemplateFile);%create presentation with the specified template

num='1';
figure
line([1 1], [100 200])
saveas(gcf,'test_plot.png')
pichandle=Picture('test_plot.png');%save figure as picture

presentationTitleSlide = add(presentationObj,'Title Slide'); %make a title slide
replace(presentationTitleSlide,'Title',['Raster Plots channel: ' num]); %rename title slide

ContentSlide = add(presentationObj,'Title and Content'); %create content slide
replace(ContentSlide,'Title','Channels XX and XX'); %replace title
replace(ContentSlide,'Content','This content bla bla') %replace content
pictureSlide = add(presentationObj,'Picture'); %Create picture slide - custom layout
for i=1:3
    figure
    line([1 1], [i i+5])
    saveas(gcf,['test_plot' num2str(i) '.png'])
    pichandle=Picture(['test_plot' num2str(i) '.png']);%save figure as picture
    replace(pictureSlide,['Picture ' num2str(i)],pichandle); %replace picture with name Picture 1
end
close(presentationObj); %close presentation to keep changes
if ispc
    winopen('testRaster.pptx'); %open presentation (WINDOWS ONLY FUNCTION)
end
