%%  MATLAB Script for Calcualting Deivation Angles and MajorAxisLengths
%  ------------------------------------------------------------------------
%  Author: Siddhartha Dhiman
%  E-Mail: sdhiman@buffalo.edu
%  Created: 07/24/2017 using Matlab R2016b
%   -----------------------------------------------------------------------

%% Clearing Workspace
clc;
clear all;
close all;
warning off;

%% Select Folder
folderCtrl = uigetdir('', 'Select Control Image Directory');
if ~isdir(folderCtrl)
    errorMessage = sprintf('Error: The following directory does not exist:\n%s', folderCtrl);
    uiwait(warndlg(errorMessage));
    return;
end

filePatternCtrl = fullfile(folderCtrl, '*.tif');
theFilesCtrl = dir(filePatternCtrl);
[upperPathCtrl, deepestFolderCtrl, ~] = fileparts(folderCtrl);

folderSch = uigetdir('', 'Select Schizophrenia Image Directory');
if ~isdir(folderSch)
    errorMessage = sprintf('Error: The following directory does not exist:\n%s', folderSch);
    uiwait(warndlg(errorMessage));
    return;
end

filePatternSch = fullfile(folderSch, '*.tif');
theFilesSch = dir(filePatternSch);
[upperPathCtrlSch, deepestFolderCtrlSch, ~] = fileparts(folderSch);

dataCtrl = getData(folderCtrl, theFilesCtrl, 1);
fprintf(1, '\n');
fprintf(1, ' Data matrix dataCtrl (Control) created \n');
fprintf(1, '\n');
fprintf(1, '\n');

dataSch = getData(folderSch, theFilesSch, 2);
fprintf(1, '\n');
fprintf(1, 'Data matrix dataSch (Schizophrenia) created \n');
fprintf(1, '\n');
fprintf(1, '\n');

%% Binning Per ROI
cVector = transpose(horzcat(dataCtrl(:).Deviation));
sVector = transpose(horzcat(dataSch(:).Deviation));
T = [cVector',sVector'];
a = min(T);
a = round(a,-1);
b = max(T);
b = round(b,-1);
space = 10;
%   Create binning edges
d = linspace(a,b,space);

% Binning loop
for i = 1:length(dataCtrl)
    binROIc = 0;
end

for i = 1:(length(d)-1)
    for k = 1:length(dataCtrl)
        clear s vc
        s = find(dataCtrl(k).Deviation >= d(i) & dataCtrl(k).Deviation ...
            < d(i+1));
        
        vc = dataCtrl(k).Deviation(s);
        
        binROIc(k,i) = length(vc);
        
    end
    for l = 1:length(dataSch)
        clear t vs
        t = find(dataSch(l).Deviation >= d(i) & dataSch(l).Deviation ...
            < d(i+1));
        vs = dataSch(l).Deviation(t);
        binROIs(l,i) = length(vs);
    end
end

%% Global Histogram
nbins = 10;
gcf = figure;
hC = histogram(cVector,nbins);
hC.FaceColor = 'r';
hC.EdgeAlpha = 0.5;
hC.EdgeColor = 'none';
hold on;
hS = histogram(sVector, nbins);
hS.FaceColor = 'b';
hS.EdgeAlpha = 0.5;
hS.EdgeColor = 'none';
hold off;
% legend('Control','Schizophrenia')
title('Histograms of Global Control and Schizophrenia Groups')
xlabel('Deviation Angle in Degrees');
ylabel('Frequency');
saveFolder = fullfile(pwd, 'Output','Global');
mkdir(saveFolder);
print(gcf, fullfile(saveFolder, sprintf('Global Histogram')),'-dpng');

%% Binning Vector and Bar
for i = 1:size(binROIc,2)
    BinVector(i,1) = sum(binROIc(:,i));
end
   
for i = 1:size(binROIs,2)
    BinVector(i,2) = sum(binROIs(:,i));
end


gcg = figure;
b = bar(BinVector);
b(1).FaceColor = [0 0.447 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
hold on
x1 = [1:1:9]-0.15;
x2 = [1:1:9]+0.15;
plot(x1,BinVector(:,1), x2,BinVector(:,2));
% legend('Control','Schizophrenia')
xlabel('Bin Number')
ylabel('Frequency')
hold off;
saveFolder = fullfile(pwd, 'Output','Global');
mkdir(saveFolder);
print(gcg, fullfile(saveFolder, sprintf('Bar Plot')),'-dpng');

%% Global CDF for Deviation
gch = figure;
g = cdfplot(cVector);
set(g,'LineWidth',2);
hold on;
h = cdfplot(sVector);
set(h,'LineWidth',2);
hold off;
xlabel('Deviation Angle in Degrees (x)');
ylabel('CDF F(x)');
% legend('Control','Schizophrenia', 'Location', 'southeast');
saveFolder = fullfile(pwd, 'Output','Global');
mkdir(saveFolder);
print(gch, fullfile(saveFolder, sprintf('CDF Plot (Deviation)')),'-dpng');

%% MajorAxisLength Binning
cVector = transpose(horzcat(dataCtrl(:).MajorAxisLength));
sVector = transpose(horzcat(dataSch(:).MajorAxisLength));
T = [cVector',sVector'];
a = min(T);
a = round(a,-1);
b = max(T);
b = round(b,-1);
space = 51;
%   Create binning edges
d = linspace(a,b,21);

% Binning loop
for i = 1:length(dataCtrl)
    binROIc = 0;
end

for i = 1:(length(d)-1)
    for k = 1:length(dataCtrl)
        clear s vc
        s = find(dataCtrl(k).MajorAxisLength >= d(i) & dataCtrl(k).MajorAxisLength ...
            < d(i+1));
        
        vc = dataCtrl(k).MajorAxisLength(s);
        
        binROIc(k,i) = length(vc);
        
    end
    for l = 1:length(dataSch)
        clear t vs
        t = find(dataSch(l).MajorAxisLength >= d(i) & dataSch(l).MajorAxisLength ...
            < d(i+1));
        vs = dataSch(l).MajorAxisLength(t);
        binROIs(l,i) = length(vs);
    end
end

%% Global CDF for MajorAxisLength
gch = figure;
g = cdfplot(cVector);
set(g,'LineWidth',2);
hold on;
h = cdfplot(sVector);
set(h,'LineWidth',2);
hold off;
xlabel('Length of Major Axis (x)');
ylabel('CDF F(x)');
% legend('Control','Schizophrenia', 'Location', 'southeast');
saveFolder = fullfile(pwd, 'Output','Global');
mkdir(saveFolder);
print(gch, fullfile(saveFolder, sprintf('CDF Plot (Length)')),'-dpng');

%% Plot No. of Events in ROIs and Distribution
for i = 1:length(dataCtrl)
    sumCtrl(i) = dataCtrl(i).TotalCells;
end

for i = 1:length(dataSch)
    sumSch(i) = dataSch(i).TotalCells;
end

xx = {'Control','Schizophrenia'};
xxC = ones(1,length(sumCtrl));
xxS = ones(1,length(sumSch)).*2;
gci = figure;
hold on;
scatter(xxC,sumCtrl,'filled','MarkerEdgeColor','k');
scatter(xxS,sumSch,'filled','MarkerEdgeColor','k');
hold off;
ylabel('No. of calretinin cells per ROI');
ax = gca;
ax.XTick = 1:numel(xx);
ax.XTickLabel = xx;
ax.XLim = [0 numel(xx)+1];
saveFolder = fullfile(pwd, 'Output','Global');
mkdir(saveFolder);
print(gci, fullfile(saveFolder, sprintf('Cell Count (per ROI)')),'-dpng');

meanCtrl = mean(sumCtrl);
stdCtrl = std(sumCtrl);
semCtrl = stdCtrl/sqrt(length(sumCtrl));
maxCtrl = 200;
minCtrl = -50;
stepCtrl = (maxCtrl-minCtrl)/100;
rangeCtrl = minCtrl:stepCtrl:maxCtrl;
pdfCtrl = normpdf(rangeCtrl, meanCtrl, stdCtrl);

meanSch = mean(sumSch);
stdSch = std(sumSch);
semSch = stdSch/sqrt(length(sumSch));
minSch = -50;
maxSch = 200;
stepSch = (maxSch-minSch)/100;
rangeSch = minSch:stepSch:maxSch;
pdfSch = normpdf(rangeSch, meanSch, stdSch);

%   Form Coordinates for Vertical Line
verticalCtrlx = [rangeCtrl(find(pdfCtrl == max(pdfCtrl)))...
    rangeCtrl(find(pdfCtrl == max(pdfCtrl)))];
verticalCtrly = [0 max(pdfCtrl)];

verticalSchx = [rangeSch(find(pdfSch == max(pdfSch))) ... 
    rangeSch(find(pdfSch == max(pdfSch)))];
verticalSchy = [0 max(pdfSch)];

gcj = figure;
hold on;
plot(minCtrl:stepCtrl:maxCtrl,pdfCtrl, 'LineWidth',2);
plot(minSch:stepSch:maxSch,pdfSch, 'LineWidth',2);
line(verticalCtrlx,verticalCtrly,'Color','b','LineStyle','--','LineWidth',2);
line(verticalSchx,verticalSchy,'Color','r','LineStyle','--','LineWidth',2);
xlabel('Number of cells per ROI');
ylabel('Probability Density');
% legend('Control','Schizophrenia');
saveFolder = fullfile(pwd, 'Output','Global');
mkdir(saveFolder);
print(gcj, fullfile(saveFolder, sprintf('Cell Count Dist')),'-dpng');



