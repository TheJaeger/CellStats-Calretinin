function data = getData(folder, theFiles, condition)
%  ------------------------------------------------------------------------
%  Author: Siddhartha Dhiman
%  E-Mail: sdhiman@buffalo.edu
%  Created: 07/24/2017 using Matlab R2016b
%  ------------------------------------------------------------------------
%  Function to generate cell counts, deviation and process length
%  into the matrix structure 'data'. This matrix contains the 
%  following information:
%       FileName:           Name of files processed
%       ROI:                ROI index number
%       TotalCells:         Total events/cells detected
%       LongProc:           Count of cells with long processes
%       ShortProc:          Count of cells with short processes
%       NoProc:             Count of cells with no processes
%       Deviation:          Row vector with devaiation of long processes
%       MajorAxisLength:    Length of all events
%   -----------------------------------------------------------------------
%   Input Arguments:
%       folder:     Path to folder containing images
%       theFiles:   Structural matrix with directory information. Can be
%                   constructed using 'dir(filePattern)', where
%                   'filePattern' is a file filtering criteria.
%       condition:  1 for Control
%                   2 for Schizophrenia
%   -----------------------------------------------------------------------
%   Output:
%       data: Structural matrix containing information useful for
%                generating biological networks.
%   -----------------------------------------------------------------------
%   Feel free to drop me an email on questions and concerns
%   -----------------------------------------------------------------------

%% Variables for Cell Detection
%   These variables can be changed to fine tune the image output
%   -----------------------------------------------------------------------
R_SE1 = 25;             % Top-Hat structuring element radius (counting)
R_SE2 = 1;              % Morphological opening structuring element radius
sq = 1;                 % Size of 'sq x sq' square matrix for marker cleaning
open_max = 100;         % Area for max. reg. clearing using 'bwareaopen'
H = 15;                 % H-maxima transform
minCellSize = 10;       % Minimum cell size in Watershed transform
sz = 25;                % Scatter plot point size

%% Structuring Elements for Cell Detection
SE1 = strel('disk', R_SE1);     % For Top-Hat (Counting)
SE2 = strel('disk', R_SE2);     % For morphological operations
SE3 = strel(ones(sq,sq));       % For marker cleaning

%% Setting up Loop Counter
cntA = 0;
cntB = 0;

%% Main Loop
for a = 1:length(theFiles)
    cntA = cntA + 1;
    baseFileName = theFiles(a).name;
    fullFileName = fullfile(folder, baseFileName);
    fprintf(1, 'Now evaluating %s. \n', baseFileName);
    imageArray = imread(fullFileName);
    [pathstr fileName ext] = fileparts(fullFileName);
    
    %% Extract Red Cahnnel
    IR = imageArray(:,:,1);
    IR(IR<31) = 0;
    
    %% Get ROI and Reference Line
    ROI = getROI(imageArray);
    line = getLine(imageArray);
    
    cntBlabel = 0;
    for b = 1:size(ROI,3)
        cntB = cntB + 1;
        cntBlabel = cntBlabel + 1;
        fprintf(1, '\ta: Working on ROI %s. \n', cntBlabel);
        
        %% Name File
        data(cntB).FileName = baseFileName;
        data(cntB).ROI = cntBlabel;
        %% Extract Red Channel
        ROImask = ROI(:,:,cntBlabel);
        maskIR = immultiply(IR,ROImask);
        
        %% Top-Hat Operation
        tophatIR = imtophat(maskIR, SE1);
        eqIR = adapthisteq(tophatIR);
        
        %% Morphological Processing
        IRo = imopen(eqIR, SE2);
        erodeIR = imerode(tophatIR, SE2);
        IRobr = imreconstruct(erodeIR, tophatIR);
        IRoc = imclose(IRo, SE2);
        IRobrd = imdilate(IRobr, SE2);
        IRobrcbr = imreconstruct(imcomplement(IRobrd), imcomplement(IRobr));
        IRobrcbr = imcomplement(IRobrcbr);
        
        %% Regional Maximas
        fgm = imextendedmax(IRobrcbr, H);
        fgm2 = imclose(fgm, SE3);
        fgm3 = imerode(fgm2, SE3);
        fgm4 = bwareaopen(fgm3, open_max);
        bw = imbinarize(IRobrcbr);
        bw1 = imclose(bw, strel('disk',5));
        
        %% Watershed Transform
        eqIR_c = imcomplement(eqIR);
        I_mod = imimposemin(eqIR_c, bw1);
        L = watershed(I_mod);
        binWatershed = L > 1;
        regsIR = regionprops(binWatershed, 'Area', 'Centroid', 'PixelIdxList',...
            'MajorAxisLength','MinorAxisLength','Orientation', 'Eccentricity');
        data(cntB).TotalCells = length(regsIR);
        
        %% Plot Ellipses with Eccentricity Thresh
        ind = find([regsIR.MajorAxisLength] > 50 & [regsIR.Eccentricity] > 0.40);
        data(cntB).LongProc = length(ind);
        data(cntB).ShortProc = length(regsIR) - length(ind);
        gcg = figure('vis', 'off');
        imshow(imageArray)
        hold on
        
        phi = linspace(0,2*pi,50);
        cosphi = cos(phi);
        sinphi = sin(phi);
        
        for k = 1:length(ind)
            xbar = regsIR(ind(k)).Centroid(1);
            ybar = regsIR(ind(k)).Centroid(2);
            
            c = regsIR(ind(k)).MajorAxisLength/2;
            d = regsIR(ind(k)).MinorAxisLength/2;
            
            theta = pi*regsIR(ind(k)).Orientation/180;
            R = [ cos(theta)   sin(theta)
                -sin(theta)   cos(theta)];
            
            xy = [c*cosphi; d*sinphi];
            xy = R*xy;
            
            x = xy(1,:) + xbar;
            y = xy(2,:) + ybar;
            
            plot(x,y,'g','LineWidth',1);
        end
        hold off
        
        if condition == 1
            saveFolder = fullfile(pwd, 'Output','Supplementary',...
                'Ellipses','Control', fileName);
            title('Detected Processes in Control');
        elseif condition == 2
            saveFolder = fullfile(pwd, 'Output','Supplementary',...
                'Ellipses','Schizophrenia', fileName);
            title('Detected Processes in Schizophrenia');
        else
            msg = 'Condition not specified or invalid when calling getData';
            error(msg)
        end
        
        mkdir(saveFolder);
        print(gcg, fullfile(saveFolder, sprintf('ROI_%d',cntBlabel)),'-dpng');
        
        %% No Processes
        indNo = find([regsIR.Eccentricity] > 0 & [regsIR.Eccentricity]...
            <= 0.30);
        data(cntB).NoProc = length(indNo);
        
        %% Convert Angles to 0 - 180 deg
        %   Line Conversion
        lineDeg = 90 - line.theta;
        
        %   Ellipse Conversion
        sFilt = regsIR(ind);
        for i = 1:length(sFilt)
            if sFilt(i).Orientation < 0;
                sFilt(i).Orientation = sFilt(i).Orientation + 180;
            else
                sFilt(i).Orientation = sFilt(i).Orientation;
            end
            sFilt(i).Deviation = abs(sFilt(i).Orientation - lineDeg);
            
            if sFilt(i).Deviation > 90;
                sFilt(i).Deviation = 180 - sFilt(i).Deviation;
            else
                sFilt(i).Deviation = sFilt(i).Deviation;
            end
        end
        
        %% Plot Histogram
        gcf = figure('vis', 'off');
        histogram([sFilt.Deviation],36);
        xlabel('Angle Deviation from Reference');
        ylabel('Frequency');
        
        if condition == 1
            saveFolder = fullfile(pwd, 'Output','Histograms','Control', fileName);
            title('Histogram of Control Angles');
        elseif condition == 2
            saveFolder = fullfile(pwd, 'Output','Histograms', 'Schizophrenia', fileName);
            title('Histogram of Schizophrenia Angles');
        else
            msg = 'Condition not specified or invalid when calling getData';
            error(msg)
        end
        
        mkdir(saveFolder);
        print(gcf, fullfile(saveFolder, sprintf('ROI_%d',cntBlabel)),'-dpng');
        
        %% Save Data
        data(cntB).ProcDetected = length(sFilt);
        for i = 1:length(sFilt)
            Deviation(i) = sFilt(i).Deviation;
        end
        data(cntB).Deviation = Deviation;
        
        for i = 1:length(regsIR)
            Length(i) = regsIR(i).MajorAxisLength;
        end
        data(cntB).MajorAxisLength = Length;
        
        %% Delete Loop Variables
        clearvars regsIR Deviation Length
    end
    clearvars ROI line IR
end
data = data;