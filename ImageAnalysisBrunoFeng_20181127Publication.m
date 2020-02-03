function [Objects] = ImageAnalysisBruno(ch1_MTR, ch2_Tom20, ch3_H, ch4_Tuj1, fileThis, FolderThisAnalysis, i)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Information|Instrument|LightSource|LightSourceType|Laser|Wavelength #1 = 561 mitotrackerRed
% Information|Instrument|LightSource|LightSourceType|Laser|Wavelength #2 = 488 Tom20
% Information|Instrument|LightSource|LightSourceType|Laser|Wavelength #3 = 405 Hoechst
% Information|Instrument|LightSource|LightSourceType|Laser|Wavelength #4 = 638 Tuj1

%vol(ch1_MTR, 0, 500)
%vol(ch2_Tom20, 0, 1000)
%vol(ch3_H, 0, 1000)
%vol(ch4_Tuj1, 0, 1000)

%% Segment mitochondria
    %% Median
    MitoDoGTom20 = imfilter(ch2_Tom20, fspecial('gaussian', 61, 1), 'symmetric') - imfilter(ch1_MTR, fspecial('gaussian', 61, 20), 'symmetric'); % vol(MitoDoGTom20, 0, 200, 'hot')
    MitoDoGthreshold = 10;
    MitoMaskTomLocal = MitoDoGTom20 > MitoDoGthreshold;%50; %vol(MitoMaskTomLocal)
    MitoDoGthresholdRefined = quantile(MitoDoGTom20(MitoMaskTomLocal), 0.5);
    MitoMaskTomLocal = MitoDoGTom20 > MitoDoGthresholdRefined;
    MitoRawthreshold = quantile(ch2_Tom20(MitoMaskTomLocal), 0.1);
    MitoMaskTomGlobal = medfilt3(ch2_Tom20) > MitoRawthreshold; % vol(MitoMaskTomGlobal)
    MitoMask = MitoMaskTomLocal & MitoMaskTomGlobal;%vol(MitoMask)
    MitoMask = bwareaopen(MitoMask, 7);%vol(MitoMask)
    
    MitoMeanMTmean = mean(ch1_MTR(MitoMask));
    
    %% Segment Nuclei
    NucleiDoG = imfilter(double(ch3_H), fspecial('gaussian', 300, 10), 'symmetric') - imfilter(double(ch3_H), fspecial('gaussian', 300, 100), 'symmetric'); %vol(NucleiDoG, 0, 200)
    NucDoGMask = NucleiDoG > 1; %vol(NucDoGMask)
    MinThreshold = 0;
    MaxThreshold = 240;
    ThresholdStep = 5;
    display = 0;
    [PreNucMask, ThresholdNuclei] = f_IterativeThresholding(ch3_H, MinThreshold, MaxThreshold, ThresholdStep, display);
    PreNucMask = PreNucMask & NucDoGMask; %vol(PreNucMask)
    ThresholdNucleiRefined = quantile(ch3_H(PreNucMask), 0.5);
    NucMask = medfilt3(ch3_H) > ThresholdNucleiRefined;
    NucMask = bwareaopen(NucMask, 20000); %vol(NucMask)
    
    %% Picnotic
    ThresholdPicnotic = 5 * ThresholdNuclei;
    PicnoticMask = medfilt3(ch3_H) > ThresholdPicnotic;
    PicnoticMask = bwareaopen(PicnoticMask, 100); % imtool(max(PicnoticMask,[],3))

    %% Collect features

    MitoLabelMatrix = bwlabeln(MitoMask, 26); %vol(MitoLabelMatrix)
    MitoObjects = regionprops('table', MitoLabelMatrix, ch1_MTR, {'Area', 'MeanIntensity'}); % doc regionprops

    % Shape
    Conn6Strel = {};
    Conn6Strel{1} = [0 0 0; 0 1 0; 0 0 0];
    Conn6Strel{2} = [0 1 0; 1 1 1; 0 1 0];
    Conn6Strel{3} = [0 0 0; 0 1 0; 0 0 0];
    Conn6Strel = logical(cat(3, Conn6Strel{:}));
    MitoErodedMask = imerode(MitoMask, Conn6Strel);
    MitoPerimMask = (MitoMask - MitoErodedMask) > 0;% vol(MitoPerimMask)

    % skeleton  
    skel = Skeleton3D(MitoMask);
    skelP = cat(3, zeros(size(skel), 'logical'), skel,zeros(size(skel), 'logical')); % vol(skelP)

    [AdjacencyMatrix, node, link] = Skel2Graph3D(skel,0); % 0 for keeping all branches  
    save([FolderThisAnalysis, filesep, num2str(i)], 'AdjacencyMatrix', 'node', 'link', 'fileThis');
    
    %% Additional feature images and vectors
    MitoBodyLabelIm = bwlabeln(MitoErodedMask, 6);

    Objects = table();
    Objects.Path = {fileThis};
    Objects.CountMito = size(MitoObjects, 1);    
    Objects.NodeDegreeVector = {sum(AdjacencyMatrix, 1)};
    Objects.MitoPixels = sum(MitoMask(:));

    % Erosion or skel derived
    Objects.MitoSkelPixels = sum(skel(:));
    Objects.MitoPerimPixels = sum(MitoPerimMask(:));
    Objects.MitoBodyPixels = sum(MitoErodedMask(:)); 
    Objects.MitoBodyCount = max(MitoBodyLabelIm(:)); % Needed for invagination feature
    Objects.MitoShapeBySurface_Norm = Objects.MitoBodyPixels / Objects.MitoPerimPixels; % Roundness feature
    Objects.MitoBodycountByMitocount_Norm = Objects.MitoBodyCount / Objects.CountMito; % Invagination feature
    
    % Nuclei features
    Objects.NucleiPixels = sum(NucMask(:));
    Objects.PicnoticPixels = sum(PicnoticMask(:));
    [~, Objects.NucleiCount] = bwlabeln(NucMask, 26);
    
    % NormalizedFeatures
    Objects.MitoSkelProportion_Norm = Objects.MitoSkelPixels / Objects.MitoPixels;
    Objects.MitoPixelsByNuc_Norm =  Objects.MitoPixels / Objects.NucleiPixels;
    Objects.MeanMitoArea_Norm = Objects.MitoPixels / Objects.CountMito;
    Objects.CountMito_Norm = Objects.CountMito / Objects.NucleiPixels;
    Objects.MitoSkelProportion_Norm = Objects.MitoSkelPixels / Objects.MitoPixels;
    Objects.MitoPerimProprtion_Norm =  Objects.MitoPerimPixels / Objects.MitoPixels;
    Objects.MitoBodyProportion_Norm = Objects.MitoBodyPixels / Objects.MitoPixels;
    Objects.MitoBodyCount_Norm = Objects.MitoBodyCount / Objects.CountMito;
    Objects.PicnosisMetric_Norm = Objects.PicnoticPixels / Objects.NucleiPixels;
    
    % Skeleton derived
    Objects.TotalNodeCount = size(node, 2);
    Objects.AverageNodePerMito_Norm = Objects.TotalNodeCount / Objects.CountMito;
    Objects.TotalLinkCount = size(link, 2);
    Objects.NodesPerMitoMean_Norm = Objects.TotalNodeCount / Objects.CountMito;
    Objects.LinksPerMitoMean_Norm = Objects.TotalLinkCount / Objects.CountMito;
    Objects.AverageNodeDegree_Norm = mean(Objects.NodeDegreeVector{:}); %average lenght of branch in pixels
    Objects.MedianNodeDegree_Norm = median(Objects.NodeDegreeVector{:});
    Objects.StdNodeDegree = std(Objects.NodeDegreeVector{:});
    Objects.MadNodeDegree = mad(Objects.NodeDegreeVector{:}, 1);
    
    % MT intensity
    Objects.MTmean = MitoMeanMTmean;
    
    %% Previews
    
    MiddlePlane = round(size(ch1_MTR, 3)/2);
    imSize = [size(ch1_MTR, 1), size(ch1_MTR, 2)];
    [BarMask, BarCenter] = f_barMask(20, 0.17195767195767195, imSize, imSize(1)-50, 50, 7); % it(BarMask)

    MitoDisplayIm3D = mat2gray(ch1_MTR); %vol(MitoDisplayIm3D)
    MitoDisplayIm = imadjust(MitoDisplayIm3D(:,:,MiddlePlane)); % imtool(MitoDisplayIm)
    PreviewMito = imoverlay(MitoDisplayIm, bwperim(MitoMask(:,:,MiddlePlane)), [1 0 0]);
    PreviewMito = imoverlay(PreviewMito, BarMask, [1 1 1]); % imtool(PreviewMito)
    
    PreviewMito3D = imoverlay(imadjust(max(MitoDisplayIm3D, [], 3)), BarMask, [1 1 1]);
    PreviewMito3D = imoverlay(PreviewMito3D, BarMask, [1 1 1]); % imtool(PreviewMito3D)
    
    PreviewMitoGamma = imoverlay(imadjust(mat2gray(max(ch2_Tom20,[],3).^2), [0 0.1], [0 1]), bwperim(max(MitoMask,[],3)), [1 0 0]);
    PreviewMitoGamma = imoverlay(PreviewMitoGamma, BarMask, [1 1 1]); % imtool(PreviewMitoGamma)
    
    PreviewMitoOnly =  imoverlay(imadjust(MitoDisplayIm3D(:,:,MiddlePlane)), BarMask, [1 1 1]); % imtool(PreviewMitoOnly)

    PreviewSkel = imoverlay(imadjust(max(MitoDisplayIm3D, [], 3)), max(skel, [], 3), [0 1 1]);
    PreviewSkel = imoverlay(PreviewSkel, BarMask, [1 1 1]); % imtool(PreviewSkel)
    
    DapiBGVec = ch3_H(~(PreNucMask));
    DapiBackground = mean(DapiBGVec);
    DapiDisplayIm3D = medfilt3(f_ImAdjust(ch3_H-DapiBackground)); % vol(DapiDisplayIm3D)

    PreviewNuclei = imoverlay(mat2gray(max(DapiDisplayIm3D, [], 3)), bwperim(max(NucMask, [], 3)), [0 0 1]);
    PreviewNuclei = imoverlay(PreviewNuclei, bwperim(max(PicnoticMask, [], 3)), [1 0 0]);
    PreviewNuclei = imoverlay( PreviewNuclei, BarMask, [1 1 1]); % imtool(PreviewNuclei)
    
    PreviewRGB = cat(3, imadjust(max(uint16(ch1_MTR), [], 3), [0 0.1], [0 1]), zeros(size(MitoDisplayIm), 'uint16'), imadjust(max(uint16(ch3_H), [], 3), [0 0.05], [0 1]));
    PreviewRGB = imoverlay(PreviewRGB, BarMask, [1 1 1]); % imtool(PreviewRGB)
    
    PreviewMitoDoG = imoverlay(imadjust(max(uint16(MitoDoGTom20),[],3)), BarMask, [1 1 1]); % imtool(PreviewMitoDoG)

    fileNameThis = regexp(fileThis, '.*\\(.*).czi', 'tokens');
    fileNameThis = fileNameThis{:}{:};
    
    filename_Mito =  [FolderThisAnalysis, filesep, 'Previews', filesep, fileNameThis, '_Mito.png'];
    imwrite(PreviewMito, filename_Mito);
    filename_Skel =  [FolderThisAnalysis, filesep, 'Previews', filesep, fileNameThis, '_Skel.png'];
    imwrite(PreviewSkel, filename_Skel);
    filename_MitoOnly =  [FolderThisAnalysis, filesep, 'Previews', filesep, fileNameThis, '_MitoOnly.png'];
    imwrite(PreviewMitoOnly, filename_MitoOnly);
    filename_Mito3D =  [FolderThisAnalysis, filesep, 'Previews', filesep, fileNameThis, '_Mito3D.png'];
    imwrite(PreviewMito3D, filename_Mito3D);
    filename_Dapi =  [FolderThisAnalysis, filesep, 'Previews', filesep, fileNameThis, '_Dapi.png'];
    imwrite(PreviewNuclei, filename_Dapi);
    filename_RGB =  [FolderThisAnalysis, filesep, 'Previews', filesep, fileNameThis, '_RGB.png'];
    imwrite(PreviewRGB, filename_RGB);
    filename_MitoDoG =  [FolderThisAnalysis, filesep, 'Previews', filesep, fileNameThis, '_MitoDoG.png'];
    imwrite(PreviewMitoDoG, filename_MitoDoG);
    filename_MitoGamma =  [FolderThisAnalysis, filesep, 'Previews', filesep, fileNameThis, '_MitoGamma.png'];
    imwrite(PreviewMitoGamma, filename_MitoGamma);

end

