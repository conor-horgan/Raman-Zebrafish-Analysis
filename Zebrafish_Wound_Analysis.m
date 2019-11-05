% -------------------------------------------------------------------------
%  Name: Zebrafish_Wound_Analysis.m
%  Version: 1.0
%  Environment: Matlab 2019a
%  Date: 22/08/2019
%  Author: Conor Horgan
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% STEP 1 - Load Raman Data
% -------------------------------------------------------------------------

% Select zebrafish Raman spectral T0 image for analysis
uiwait(msgbox('Select the zebrafish T0 Raman spectral data to import for analysis', 'Select Spectral Data', 'modal'));
Zebrafish_Image_Input_T0 = uigetvariables({'Please select input zebrafish T0 spectral image'});
Zebrafish_Image_T0 = Zebrafish_Image_Input_T0{1};

% Select zebrafish Raman spectral T1 image for analysis
uiwait(msgbox('Select the zebrafish T1 Raman spectral data to import for analysis', 'Select Spectral Data', 'modal'));
Zebrafish_Image_Input_T1 = uigetvariables({'Please select input zebrafish T1 spectral image'});
Zebrafish_Image_T1 = Zebrafish_Image_Input_T1{1};

% Select zebrafish Raman spectral T2 image for analysis
uiwait(msgbox('Select the zebrafish T2 Raman spectral data to import for analysis', 'Select Spectral Data', 'modal'));
Zebrafish_Image_Input_T2 = uigetvariables({'Please select input zebrafish T2 spectral image'});
Zebrafish_Image_T2 = Zebrafish_Image_Input_T2{1};

% Select spectral axis for plotting
uiwait(msgbox('Select the spectral axis for plotting', 'Select Spectral Axis', 'modal'));
X_Axis_Input = uigetvariables({'Please select spectral axis'});
X_Axis = X_Axis_Input{1};


% -------------------------------------------------------------------------
% STEP 3 - Combine input data into a single matrix for analysis
% -------------------------------------------------------------------------
Data = Merge_Images(Zebrafish_Image_T0, Zebrafish_Image_T1, Zebrafish_Image_T2);


% -------------------------------------------------------------------------
% STEP 5 - Normalize input data
% -------------------------------------------------------------------------
Data_Normalised = normaliz(Data,0,1);

Zebrafish_Image_T0_norm = Data_Normalised(1:(size(Zebrafish_Image_T0,1) * size(Zebrafish_Image_T0,2)),:);
Zebrafish_Image_T1_norm = Data_Normalised((size(Zebrafish_Image_T0,1) * size(Zebrafish_Image_T0,2)) + 1:(size(Zebrafish_Image_T0,1) * size(Zebrafish_Image_T0,2)) * 2,:);
Zebrafish_Image_T2_norm = Data_Normalised((size(Zebrafish_Image_T0,1) * size(Zebrafish_Image_T0,2) * 2) + 1:(size(Zebrafish_Image_T0,1) * size(Zebrafish_Image_T0,2)) * 3,:);
Zebrafish_Image_T0_norm = reshape(Zebrafish_Image_T0_norm, size(Zebrafish_Image_T0,1), size(Zebrafish_Image_T0,2), size(Zebrafish_Image_T0,3));
Zebrafish_Image_T1_norm = reshape(Zebrafish_Image_T1_norm, size(Zebrafish_Image_T1,1), size(Zebrafish_Image_T1,2), size(Zebrafish_Image_T1,3));
Zebrafish_Image_T2_norm = reshape(Zebrafish_Image_T2_norm, size(Zebrafish_Image_T2,1), size(Zebrafish_Image_T2,2), size(Zebrafish_Image_T2,3));


% -------------------------------------------------------------------------
% STEP 6 - Remove spectra with cosmic rays from analysis
% -------------------------------------------------------------------------
Threshold = 0.017;

Data_Normalised_Cosmic = Data_Normalised;
Data_Diff = diff(diff(Data_Normalised));
tmp = find(max(Data_Diff(:,:)') > Threshold);
Data_Normalised_Cosmic([tmp], :)=[];
Data_Normalised_Cosmic = sgolayfilt(Data_Normalised_Cosmic, 1, 7);


% -------------------------------------------------------------------------
% STEP 7 - Perform VCA on normalized cosmic-less spectra
% -------------------------------------------------------------------------
% Perform VCA with 6 endmembers
% VCA: J. Nascimento and J. Dias, "Vertex Component Analysis: A fast algorithm to unmix hyperspectral data", 
% IEEE Transactions on Geoscience and Remote Sensing, vol. 43, no. 4, pp. 898-910, 2005.
% Code available from: http://www.lx.it.pt/~bioucas/code.htm
[Ae, indice, Rp] = vca2(Data_Normalised_Cosmic', 'Endmembers', 6);

% View VCA Endmembers
Display_Multiple_Spectra(X_Axis, Ae(:,1), Ae(:,2), Ae(:,3), Ae(:,4), Ae(:,5), Ae(:,6));

% Select VCA endmembers for regression
uiwait(msgbox('Select VCA endmembers for regression', 'Select VCA Endmembers', 'modal'));
Prompt = {'Enter endmember index 1:','Enter endmember index 2:', 'Enter endmember index 3:', 'Enter endmember index 4:'};
Dlgtitle = 'Endmember Indices';
Dims = [1 50];
Definput = {'','','',''};
Endmember_Indices = inputdlg(Prompt,Dlgtitle,Dims,Definput);


% -------------------------------------------------------------------------
% STEP 8 - Perform least-squares regression with selected VCA endmembers
% -------------------------------------------------------------------------
for i = 1:size(Data_Normalised,1)
    Data_Normalised_lsq(i,:) = lsqnonneg(Ae(:, [str2num(Endmember_Indices{1}) str2num(Endmember_Indices{2}) str2num(Endmember_Indices{3}) str2num(Endmember_Indices{4})]), Data_Normalised(i,:)');
end

% -------------------------------------------------------------------------
% STEP 9 - Unmerge image data into original matrices
% -------------------------------------------------------------------------
[Zebrafish_Image_T0_Output_1, Zebrafish_Image_T1_Output_1, Zebrafish_Image_T2_Output_1] = Unmerge_Images(Data_Normalised_lsq,1,Zebrafish_Image_T0, Zebrafish_Image_T1, Zebrafish_Image_T2);
[Zebrafish_Image_T0_Output_2, Zebrafish_Image_T1_Output_2, Zebrafish_Image_T2_Output_2] = Unmerge_Images(Data_Normalised_lsq,2,Zebrafish_Image_T0, Zebrafish_Image_T1, Zebrafish_Image_T2);
[Zebrafish_Image_T0_Output_3, Zebrafish_Image_T1_Output_3, Zebrafish_Image_T2_Output_3] = Unmerge_Images(Data_Normalised_lsq,3,Zebrafish_Image_T0, Zebrafish_Image_T1, Zebrafish_Image_T2);
[Zebrafish_Image_T0_Output_4, Zebrafish_Image_T1_Output_4, Zebrafish_Image_T2_Output_4] = Unmerge_Images(Data_Normalised_lsq,4,Zebrafish_Image_T0, Zebrafish_Image_T1, Zebrafish_Image_T2);


% -------------------------------------------------------------------------
% STEP 10 - View heatmap images for each endmember
% -------------------------------------------------------------------------
Display_Multiple_Images(Zebrafish_Image_T0_Output_1, Zebrafish_Image_T1_Output_1, Zebrafish_Image_T2_Output_1);
Display_Multiple_Images(Zebrafish_Image_T0_Output_2, Zebrafish_Image_T1_Output_2, Zebrafish_Image_T2_Output_2);
Display_Multiple_Images(Zebrafish_Image_T0_Output_3, Zebrafish_Image_T1_Output_3, Zebrafish_Image_T2_Output_3);
Display_Multiple_Images(Zebrafish_Image_T0_Output_4, Zebrafish_Image_T1_Output_4, Zebrafish_Image_T2_Output_4);


% -------------------------------------------------------------------------
% STEP 11 - View heatmap images for final image timepoint
% -------------------------------------------------------------------------
Display_Multiple_Images(Zebrafish_Image_T2_Output_1, Zebrafish_Image_T2_Output_2, Zebrafish_Image_T2_Output_3, Zebrafish_Image_T2_Output_4);


% -------------------------------------------------------------------------
% STEP 12 - Select VCA endmember corresponding to wounded tissue for analysis
% -------------------------------------------------------------------------
uiwait(msgbox('Select VCA endmember corresponding to wounded tissue for analysis', 'Select Wounded Tissue VCA Endmember', 'modal'));
Wound_VCA_Input = uigetvariables({'Please select wounded tissue VCA endmember, Zebrafish_Image_T2_Output_X'});
Wound_VCA = Wound_VCA_Input{1};


% -------------------------------------------------------------------------
% STEP 13 - Extract wounded and non-wounded regions
% -------------------------------------------------------------------------
Crop_Start = 13;
Fish = Zebrafish_Image_T2(Crop_Start:end,:,:);

Wound_Component_Peak_Threshold = im2bw(Wound_VCA, 0.00032);
Non_Wound_Component_Peak_Threshold = imcomplement(Wound_Component_Peak_Threshold);

figure(7)
subplot(1,2,1)
imshow(Wound_Component_Peak_Threshold)
subplot(1,2,2)
imshow(Non_Wound_Component_Peak_Threshold)

[Wound_Component_X,Wound_Component_Y] = find(Wound_Component_Peak_Threshold > 0);
Wound_Coordinates = [Wound_Component_X,Wound_Component_Y];
Wound_Coordinates_Crop = Wound_Coordinates(:,1) < Crop_Start;
Wound_Coordinates_New = Wound_Coordinates;
Wound_Coordinates_New(Wound_Coordinates_Crop,:) = [];

[Non_Wound_Component_X,Non_Wound_Component_Y] = find(Non_Wound_Component_Peak_Threshold > 0);
Non_Wound_Coordinates = [Non_Wound_Component_X,Non_Wound_Component_Y];
Non_Wound_Coordinates_Crop = Non_Wound_Coordinates(:,1) < Crop_Start;
Non_Wound_Coordinates_New = Non_Wound_Coordinates;
Non_Wound_Coordinates_New(Non_Wound_Coordinates_Crop,:) = [];


% -------------------------------------------------------------------------
% STEP 14 - Extract wounded and non-wounded spectra
% -------------------------------------------------------------------------
Wound_Component_Spectra = zeros(length(Wound_Coordinates_New),length(X_Axis));

for j = 1:length(Wound_Coordinates_New)
    Wound_Component_Spectra(j,:) = Zebrafish_Image_T2_norm(Wound_Coordinates_New(j,1), Wound_Coordinates_New(j,2), :);
end

Non_Wound_Component_Spectra = zeros(length(Non_Wound_Coordinates_New),length(X_Axis));

for j = 1:length(Non_Wound_Coordinates_New)
    Non_Wound_Component_Spectra(j,:) = Zebrafish_Image_T2_norm(Non_Wound_Coordinates_New(j,1), Non_Wound_Coordinates_New(j,2), :);
end


% -------------------------------------------------------------------------
% STEP 15 - Remove cosmic peak spectra from wounded and non-wounded spectra
% -------------------------------------------------------------------------
Threshold_2 = 0.015;

Wound_Component_Spectra_Cosmic = Wound_Component_Spectra;
Data_Diff_Wound = diff(diff(Wound_Component_Spectra_Cosmic));
%plot(Data_Diff_Wound)
tmp2 = find(max(Data_Diff_Wound(:,:)')>Threshold_2);
Wound_Component_Spectra_Cosmic([tmp2],:)=[];
Wound_Component_Spectra_cosmic = Wound_Component_Spectra_Cosmic;

Non_Wound_Component_Spectra_Cosmic = Non_Wound_Component_Spectra;
Data_Diff_Non_Wound = diff(diff(Non_Wound_Component_Spectra_Cosmic));
%plot(Data_Diff_Non_Wound)
tmp3 = find(max(Data_Diff_Non_Wound(:,:)')>Threshold_2);
Non_Wound_Component_Spectra_Cosmic([tmp3],:)=[];

% -------------------------------------------------------------------------
% STEP 16 - Calculate mean and standard deviation wounded and non-wounded
% spectra and the differnce spectrum between wounded and non-wounded 
% -------------------------------------------------------------------------
Mean_Non_Wound_Spectrum = sgolayfilt(mean(Non_Wound_Component_Spectra_Cosmic),1,7);
Mean_Wound_Spectrum = sgolayfilt(mean(Wound_Component_Spectra_Cosmic),1,7);
STD_Non_Wound_Spectrum = std(Non_Wound_Component_Spectra_Cosmic);
STD_Wound_Spectrum = std(Wound_Component_Spectra_Cosmic);

Difference_Spectrum = Mean_Non_Wound_Spectrum - Mean_Wound_Spectrum;

Plot_Multiple_Raman_Spectrum(X_Axis,Mean_Non_Wound_Spectrum, Mean_Wound_Spectrum)
