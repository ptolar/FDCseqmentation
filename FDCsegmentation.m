
% main script to segment FDC networks in 3d LN images
% subsegments FDC networks into concentric rings to quantify antigen localisation
% requires MIJ to display images (http://bigwww.epfl.ch/sage/soft/mij/)

%before running the script, load CD35 image stack as a matrix (x,y,z) named imageMarker and antigen-IC image stack as imageSignal

%% set up
shrink = 4; % reduction factor for speed
cellSize = 20; %40
follicleSize = 400; %200
th = 90; % threshold for follicle detection 5 if 8-b and 90 if 12b
splitControl = 1; % 0 more splitting to 2 less splitting
numRings = 6;
fdcTh = 2; % threshold for FDC segmentation: how many times  over background

if ~ exist('imageSignal', 'var')
    imageSignal = zeros(size(imageMarker));
end

imageMarker = padarray(imageMarker, [0,0,1]);
imageSignal = padarray(imageSignal, [0,0,1]);

%% subsample image?
if size(imageMarker,1)>1000
    imageMarkerS = zeros(ceil(size(imageMarker,1)/shrink), ceil(size(imageMarker,2)/shrink), size(imageMarker,3));
    imageSignalS = zeros(ceil(size(imageSignal,1)/shrink), ceil(size(imageSignal,2)/shrink), size(imageSignal,3));
    for m = 1:size(imageMarker,3)
        imageMarkerS(:,:,m) = imresize(imageMarker(:,:,m), 1/shrink);
        imageSignalS(:,:,m) = imresize(imageSignal(:,:,m), 1/shrink);
    end
    imageMarker = imageMarkerS;
    imageSignal = imageSignalS;
    clear imageMarkerS
    clear imageSignalS
    
end

cellSize = round(cellSize/shrink);
follicleSize = round(follicleSize/shrink);
%% filter
imageMarker = padarray(imageMarker, [follicleSize, follicleSize, 0]);
mask = zeros(size(imageMarker));
for m = 1: size(imageMarker,3)
    mask(:,:,m) = bpass(imageMarker(:,:,m), cellSize, follicleSize);
end
mask = mask(follicleSize+1:end-follicleSize, follicleSize+1:end-follicleSize, :);
imageMarker = imageMarker(follicleSize+1:end-follicleSize, follicleSize+1:end-follicleSize, :);

MIJ.createImage('mask', mask, 1);
MIJ.createImage('imageMarker', imageMarker, 1);

%% segment follicles
follicleMask = mask>th;

% some watershed splitting
d = -bwdist(~follicleMask);
dmask = imextendedmin(d,splitControl);
d = imimposemin(d, dmask);
w = watershed(d);
follicleMask(w==0) = 0;

%% remove follicles touching image borders and small ones
% follicleMask = imclearborder(follicleMask);
% props = regionprops(follicleMask, 'Area');
% volume = [props(:).Area];
% ind = 10 * volume > pi*(size(imageMarker,3)/2)*(follicleSize/2)^2;
% ind = find(~ind);
% follicleMask = bwlabeln(follicleMask);
% for m=1:numel(ind)
%     follicleMask(follicleMask==ind(m)) = 0;
% end
% follicleMask = follicleMask>0;

%% label follicles
follicleMask = bwlabeln(follicleMask);
numFollicles = max(follicleMask(:));

%% subsegment each follicle and get data from rings
dataSignal = nan(numFollicles, numRings);
dataSignalOnMarker = nan(numFollicles, numRings);
dataMarker = nan(numFollicles, numRings);
dataCoverage = nan(numFollicles,1);
for m = 1:numFollicles
    fol = follicleMask == m;
    rings = centripetalSegmentation(fol, numRings);
    if ~isempty(rings)
        for r = 1:numRings
            dataSignal(m,r) = mean(imageSignal(rings==r));
            dataMarker(m,r) = mean(imageMarker(rings==r));
            dataSignalOnMarker(m,r) = mean(imageSignal(rings==r & imageMarker>fdcTh*dataMarker(m,1)));
        end
        dataCoverage(m) = sum(imageSignal(fol & imageSignal>fdcTh*dataSignal(m,1) & imageMarker>fdcTh*dataMarker(m,1))) ...
            ./ sum(imageMarker(fol & imageMarker>dataMarker(m,1)));
    end
end

%% put results into tables
varNames = cell(1, 2*numRings+1);
varNames{1} = 'Follicle';
for m = 1:numRings
    varNames{m+1} = ['Ring_', num2str(m), '_Marker_Intensity'];
end
for m = 1:numRings
    varNames{numRings+m+1} = ['Ring_', num2str(m), '_Signal_Intensity'];
end

data = [[1:numFollicles]', dataMarker, dataSignal];
data = array2table(data, 'VariableNames', varNames);

bckgMarker = repmat(dataMarker(:,1), [1,numRings]);
bckgSignal = repmat(dataSignal(:,1), [1,numRings]);
meanMarker = repmat(mean(dataMarker-bckgMarker,2), [1,numRings]);
meanSignal = repmat(mean(dataSignal-bckgSignal,2), [1,numRings]);
meanSignalOnMarker = repmat(mean(dataSignalOnMarker-bckgSignal,2), [1,numRings]);
dataMarkerNorm = (dataMarker-bckgMarker)./meanMarker;
dataSignalNorm = (dataSignal-bckgSignal)./meanSignal;
dataSignalOnMarkerNorm = (dataSignalOnMarker-bckgSignal)./meanSignalOnMarker;

dataNorm = [[1:numFollicles]', dataMarkerNorm, dataSignalNorm];
dataNorm = array2table(dataNorm, 'VariableNames', varNames);

% dataRatio was used as final data
dataRatio = dataSignalNorm(:,2:end)./dataMarkerNorm(:,2:end);
dataRatio = array2table(dataRatio, 'VariableNames', varNames(2:numRings));



%% show images
MIJ.createImage('imageMarker', imageMarker, 1);
MIJ.createImage('follicleMask', follicleMask, 1);
MIJ.run('Merge Channels...', 'c5=imageMarker c6=follicleMask create');


