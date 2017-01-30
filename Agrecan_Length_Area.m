%Please read more about this code and how to use it in README.txt file

%this script should analyse the particles length and area.
close all

%please select image name and a path to save the outputs here.

path = 'analysed_images_4\';


Miji();
MIJ.run('Open...');
%extracts the image.
I = MIJ.getCurrentImage;
E = uint8(I);

%Perform the processing in ImageJ
MIJ.run('Run...'); %Here it will promt user to apply a macro to image.
I2 = MIJ.getCurrentImage;
imageas_full=char(MIJ.getCurrentTitle);
C = strsplit(imageas_full,'.');
imageas =C{1,1};
d=MIJ.version;

E2 = uint8(I2);
originalImage=E2;


binaryImage = im2bw(originalImage);
imagesc(binaryImage);
r=int16(binaryImage);



labeledImage = bwlabel(binaryImage, 8);
coloredLabels = label2rgb (labeledImage, 'hsv', 'k', 'shuffle');
blobMeasurements = regionprops(labeledImage, originalImage, 'all');


%create an Image with a blob number anotations
t=I.*r;
fig=figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(t);
hold on

for k = 1 : length(blobMeasurements)         
	thisCentroid = [blobMeasurements(k).Centroid(1), blobMeasurements(k).Centroid(2)];
	message = num2str((k));
    %add a number to each Agrecan molecules for further validation
	text(thisCentroid(1), thisCentroid(2), message, 'Color', 'y');
	
end

%%%%%%%%%%%%%%%
%calculates the length of the protein
skeleton = bwmorph(binaryImage,'thin',Inf);
edtImage = bwdist(~binaryImage);
distanceFromEdge = edtImage .* single(skeleton);

%select the area of each protein and then calculate the length of the
%protein - skeletal distanceFromEdge image pixel number greater than 0 + the 
%average distance number. http://www.mathworks.com/matlabcentral/answers/43506-polygon-width-and-centerline
%http://uk.mathworks.com/help/images/ref/bwmorph.html



%This perfors an analysis of each of Agrecan molecules
ProteinAnalysis=[];

for r1=1:1:length(blobMeasurements)

Area=blobMeasurements(r1).Area;
    
Pixels=blobMeasurements(r1).PixelIdxList;

SkeletPix=distanceFromEdge(Pixels);

%remove the part where there is no line visible
SingleLinePix=find(SkeletPix>0);

%finds the length of the line in pixels

ProtInLength=length(SingleLinePix);

%fixes the error of the length adding the average length of closest point
%distance.

%length error correction

Error= round(mean(SkeletPix(SingleLinePix)));
ProtLength=ProtInLength+(2*Error);

%puts everything in one matrix 

ProteinAnalysis(r1,1)=r1;
ProteinAnalysis(r1,2)=Area;
ProteinAnalysis(r1,3)=ProtLength;



end


%Skeleton_Image_display

t=I.*r;
fig2=figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(skeleton)
hold on


for k = 1 : length(blobMeasurements)          % Loop through all blobs adding a number ID of agrecan.
	thisCentroid = [blobMeasurements(k).Centroid(1), blobMeasurements(k).Centroid(2)];
	message = num2str((k));
	text(thisCentroid(1), thisCentroid(2), message, 'Color', 'y');
	
end

%%%%%%%%%%%%%%%
%calculates the length of the protein
skeleton = bwmorph(binaryImage,'thin',Inf);
edtImage = bwdist(~binaryImage);
distanceFromEdge = edtImage .* single(skeleton);

%Prepere Image for an export

myStyle = hgexport('factorystyle');
myStyle.Resolution = 300;

%Here the analysed results are exported
mkdir(strcat(path,imageas));
hgexport(fig2, strcat(path,imageas,'\ImageSkeleton'), myStyle, 'Format', 'jpeg');
hgexport(fig, strcat(path,imageas,'\Image'), myStyle, 'Format', 'jpeg');
dlmwrite(strcat(path,imageas,'\ProteinAreaLength.txt'),ProteinAnalysis,'delimiter', '\t');



%code written by Matiss Ozols 30/01/2017