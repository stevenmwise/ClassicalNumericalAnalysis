usamap('Tennessee')

tenni = shaperead('usastatehi', 'UseGeoCoords', true,...
  'Selector',{@(name) strcmpi(name,'Tennessee'), 'Name'});

geoshow(tenni, 'FaceColor', [0.9290 0.6940 0.1250])

textm(tenni.LabelLat, tenni.LabelLon, tenni.Name,...
  'HorizontalAlignment', 'center')

lenLon = length(tenni.Lon);
lenLat = length(tenni.Lat);

close all;
initmsh();

node = zeros(lenLon-1,2);

node(:,1) = tenni.Lon(1:lenLon-1);
node(:,2) = tenni.Lat(1:lenLat-1);

edge = zeros(lenLon-1,2);

edge(:,1) = 1:lenLon-1;
edge(:,2) = 2:lenLon;
edge(lenLon-1,2) = 1;

opts.kind = 'delaunay';
opts.rho2 = +1.0 ;

[vert,etri,tria,tnum] = refine2(node,edge,[],opts);

hf = figure(1);
patch('faces',tria(:,1:3),'vertices',vert, ...
   'facecolor','w', ...
   'edgecolor',[.2,.2,.2]) ;
 
hold on; axis image off;
patch('faces',edge(:,1:2),'vertices',node, ...
  'facecolor','w', ...
  'edgecolor',[.1,.1,.1], ...
  'linewidth',1.0);

exportgraphics(hf,'tennessee.eps')
exportgraphics(hf,'tennessee.pdf')

