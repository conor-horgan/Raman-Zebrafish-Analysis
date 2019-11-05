% -------------------------------------------------------------------------
%  Name: Zebrafish_Granuloma_Analysis.m
%  Version: 1.0
%  Environment: Matlab 2019a
%  Date: 15/08/2019
%  Author: Conor Horgan
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
%  1. Select variables for analysis
% -------------------------------------------------------------------------
% Select zebrafish Raman spectral image for analysis
Zebrafish_Image_Input = uigetvariables({'Please select input zebrafish spectral image'});
Zebrafish_Image = Zebrafish_Image_Input{1};

% Select spectral index range for analysis
Prompt = {'Enter start index range:','Enter end index range:'};
Dlgtitle = 'Spectral Index Range';
Dims = [1 50];
Definput = {'964','974'};
Spectral_Range = inputdlg(Prompt,Dlgtitle,Dims,Definput);


% -------------------------------------------------------------------------
%  2. Set paramter values
% -------------------------------------------------------------------------
C = gray(128);
C_size = size(C,1);

Peak_Threshold = 0.0019;
BW_Threshold = 0.3;
BW_PixelSize = 5;


% -------------------------------------------------------------------------
%  3. Calculate Peak Intensity Image
% -------------------------------------------------------------------------
Peak = double(Zebrafish_Image(:,:,str2num(Spectral_Range{1}):str2num(Spectral_Range{2})));
Peak_Sum = sum(Peak,3);
Peak_Scaled = round(interp1(linspace(min(Peak_Sum(:)),max(Peak_Sum(:)),C_size),1:C_size,Peak_Sum));
Peak_Image = interp2(1:size(Peak_Scaled,2), 1:size(Peak_Scaled,1), Peak_Scaled, 1:0.05:size(Peak_Scaled,2), (1:0.05:size(Peak_Scaled,1))');
% imwrite(mat2gray(Peak_Image), 'Granuloma_Peak.png');


% -------------------------------------------------------------------------
%  3. Generate 'Whole Zebrafish and Granuloma' Image
% -------------------------------------------------------------------------
Whole_Fish = double(Zebrafish_Image(:,:,943:1054));
Whole_Fish_Sum = sum(Whole_Fish,3);
Whole_Fish_Scaled = round(interp1(linspace(min(Whole_Fish_Sum(:)),max(Whole_Fish_Sum(:)),C_size),1:C_size,Whole_Fish_Sum));
Whole_Fish_Image = interp2(1:size(Whole_Fish_Scaled,2), 1:size(Whole_Fish_Scaled,1), Whole_Fish_Scaled, 1:0.05:size(Whole_Fish_Scaled,2), (1:0.05:size(Whole_Fish_Scaled,1))');
% imwrite(mat2gray(Whole_Fish_Image), 'Whole_Fish_.png');


% -------------------------------------------------------------------------
%  4. Apply Peak Intensity Threshold
% -------------------------------------------------------------------------
Peak_Threshold = Peak_Sum;
for j = 1:size(Peak_Threshold,1)
    for k = 1:size(Peak_Threshold,2)
        if Peak_Threshold(j,k) < Peak_Threshold
            Peak_Threshold(j,k) = 0;
        end
    end
end


% -------------------------------------------------------------------------
%  5. Extract Peak Thresholded Pixel Spectra
% -------------------------------------------------------------------------
[Peak_X,Peak_Y] = find(Peak_Threshold > 0);

Peak_Spectra = zeros(length(Peak_X),1383);

for j = 1:length(Peak_X)
    Peak_Spectra(j,:) = Zebrafish_Image(Peak_X(j), Peak_Y(j), :);
end


% -------------------------------------------------------------------------
%  6. Isolate Thresholded Regions for Binary Image
% -------------------------------------------------------------------------
Threshold_Image = mat2gray(Peak_Threshold);
BW_Image = im2bw(Threshold_Image, BW_Threshold);
BW_Image_Crop = bwareaopen(BW_Image,BW_PixelSize);
stats = regionprops(BW_Image_Crop,'Centroid','Area', 'PixelList');


% -------------------------------------------------------------------------
%  7. Assign Thresholded Spectra to Thresholded Image Regions
% -------------------------------------------------------------------------
Centroid_Membership = zeros(size(Peak_Spectra,1),1);

for i = 1:size(Peak_Spectra,1)
    Coordinate = [Peak_Y(i) Peak_X(i)];
    for j = 1:size(stats,1)
        if ismember(Coordinate, stats(j).PixelList,'rows')
            Centroid_Membership(i) = j;
        end
    end
end

Centroid_Membership = nonzeros(Centroid_Membership);

Centroid_Spectra = cell(size(stats,1),1);
for j = 1:size(stats,1)
    Centroid_Spectra{j,1} = Peak_Spectra(find(Centroid_Membership == j),:);
end


% -------------------------------------------------------------------------
%  8. Calculate Mean Spectra for each Thresholded Image Region
% -------------------------------------------------------------------------
Mean_Centroid_Spectra = zeros(size(stats,1),1383);
for j = 1:size(stats,1)
    if size(Centroid_Spectra{j,1},1) > 1
        Mean_Centroid_Spectra(j,:) = mean(Centroid_Spectra{j,1});
    else
        Mean_Centroid_Spectra(j,:) = Centroid_Spectra{j,1};
    end
end

All_Centroid_Spectra = [];
for j = 1:size(stats,1)
    All_Centroid_Spectra = [All_Centroid_Spectra; Centroid_Spectra{j,1}];
end

[Centroid_Membership_sorted, Centroid_Membership_order] = sort(Centroid_Membership);
All_Centroid_Spectra_sorted = All_Centroid_Spectra(Centroid_Membership_order,:);


% -------------------------------------------------------------------------
%  9. Define Colormaps for Labelled Image
% -------------------------------------------------------------------------
Colormap = [240 0 130; 250 60 60; 240 130 40; 230 175 45; 230 220 50; 160 230 50; 0 220 0; 0 210 140; 0 200 200; 0 160 225; 30 60 255; 110 0 220; 160 0 200; 255 0 0]/255;
Colormap_LightBlue = [183 228 255; 183 228 255; 183 228 255; 183 228 255; 183 228 255; 183 228 255; 183 228 255; 183 228 255; 183 228 255; 183 228 255; 183 228 255; 183 228 255; 183 228 255; 183 228 255; 183 228 255]/255;
Colormap_MedBlue = [79 189 255; 79 189 255; 79 189 255; 79 189 255; 79 189 255; 79 189 255; 79 189 255; 79 189 255; 79 189 255; 79 189 255; 79 189 255; 79 189 255; 79 189 255; 79 189 255; 79 189 255; 79 189 255]/255;
Colormap_DarkBlue = [0 72 115; 0 72 115; 0 72 115; 0 72 115; 0 72 115; 0 72 115; 0 72 115; 0 72 115; 0 72 115; 0 72 115; 0 72 115; 0 72 115; 0 72 115; 0 72 115; 0 72 115; 0 72 115]/255;

Colormap_LightRed = [255 183 183; 255 183 183; 255 183 183; 255 183 183; 255 183 183; 255 183 183; 255 183 183; 255 183 183; 255 183 183; 255 183 183; 255 183 183; 255 183 183; 255 183 183; 255 183 183; 255 183 183]/255;
Colormap_MedRed = [255 2 2; 255 2 2; 255 2 2; 255 2 2; 255 2 2; 255 2 2; 255 2 2; 255 2 2; 255 2 2; 255 2 2; 255 2 2; 255 2 2; 255 2 2; 255 2 2; 255 2 2; 255 2 2]/255;
Colormap_DarkRed = [115 0 0; 115 0 0; 115 0 0; 115 0 0; 115 0 0; 115 0 0; 115 0 0; 115 0 0; 115 0 0; 115 0 0; 115 0 0; 115 0 0; 115 0 0; 115 0 0; 115 0 0; 115 0 0]/255;

[labeledImage numberOfBlobs] = bwlabel(BW_Image_Crop, 8);
coloredLabelsImage = label2rgb(labeledImage, Colormap, 'k');


% -------------------------------------------------------------------------
%  10. Display Image Processing Results
% -------------------------------------------------------------------------
figure(1)
suptitle('Zebrafish Granuloma Analysis')
subplot(2,2,1)
imshow(Threshold_Image)
title('Thresholded Image')

subplot(2,2,2)
imshow(BW_Image)
title('Binary Thresholded Image')

subplot(2,2,3)
imshow(coloredLabelsImage)
title('Colored Labelled Image')

subplot(2,2,4)
set(groot, 'defaultAxesColorOrder', Colormap)
plot(1:1383, Mean_Centroid_Spectra)
ylim([0 1.1*max(max(Mean_Centroid_Spectra))])
xlabel('Spectral Axis Index')
ylabel('Intensity (a.u.)')
title('Mean Centroid Spectra')
set(groot,'defaultAxesColorOrder','remove')

