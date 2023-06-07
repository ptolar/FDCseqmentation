function rings = centripetalSegmentation(object, segments)
% segments object mask into circular centripetal rings
% segments is the number of rings to use


if nargin<2
    segments = 4;
end

% calculate distance from perimeter for each pixel
if ndims(object) == 2
    % 2D case
    edge = bwperim(object);
else
    % 3D case
    edge = object - imerode(object, true(3,3,3));
end
dist = bwdist(edge);

% bin distance values
binEdges = 0 : max(dist(object))/segments : max(dist(object));
if ~isempty(binEdges)
    binEdges(end) = binEdges(end)+1;
    [~, rings] = histc(dist, binEdges);
    rings(~object) = 0;
else
    rings = [];
end



