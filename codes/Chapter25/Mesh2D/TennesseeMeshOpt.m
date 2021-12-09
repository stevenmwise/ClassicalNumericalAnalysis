usamap('Tennessee')

tenni = shaperead('usastatehi', 'UseGeoCoords', true,...
  'Selector',{@(name) strcmpi(name,'Tennessee'), 'Name'});

%geoshow(tenni, 'FaceColor', [0.9290 0.6940 0.1250])

lenLon = length(tenni.Lon);
lenLat = length(tenni.Lat);

close all;
initmsh( );

node = zeros(lenLon-2,2);

node(1:lenLon-2,1) = tenni.Lon(1:lenLon-2);
node(1:lenLon-2,2) = tenni.Lat(1:lenLat-2);

edge = zeros(lenLon-2,2);

edge(1:lenLon-2,1) = 1:lenLon-2;
edge(1:lenLon-2,2) = 2:lenLon-1;
edge(  lenLon-2,2) = 1;

olfs.dhdx = +0.15;

[vlfs,tlfs,hlfs] = lfshfn2(node,edge,[ ],olfs);

[slfs] = idxtri2(vlfs,tlfs);

% Generate the mesh:
    
hfun = @trihfn2;

[vert,etri,tria,tnum] = refine2(node,edge,[ ],[ ],hfun, ...
  vlfs,tlfs,slfs,hlfs);

hf1 = figure(1);
    
patch('faces',tria(:,1:3),'vertices',vert, ...
  'facecolor','w', ...
  'edgecolor',[.2,.2,.2]);
    
hold on; 
axis image off;

patch('faces',edge(:,1:2),'vertices',node, ...
  'facecolor','w', ...
  'edgecolor',[.1,.1,.1], ...
  'linewidth',0.1);

% Optimize the mesh:

[vnew,enew,tnew,tnum] = smooth2(vert,etri,tria,tnum);

hf2 = figure(2);
patch('faces',tnew(:,1:3),'vertices',vnew, ...
  'facecolor',[0.9290 0.6940 0.1250], ...
  'edgecolor',[.2,.2,.2]) ;

hold on;
axis image off;

patch('faces',edge(:,1:2),'vertices',node, ...
  'facecolor',[0.9290 0.6940 0.1250], ...
  'edgecolor',[.1,.1,.1], ...
  'linewidth',0.1);

exportgraphics(hf2,'tennessee.eps','ContentType','vector')
exportgraphics(hf2,'tennessee.pdf','ContentType','vector')

%hvrt = trihfn2(vert,vlfs,tlfs,slfs,hlfs) ;
%hnew = trihfn2(vnew,vlfs,tlfs,slfs,hlfs) ;

%tricost(vert,etri,tria,tnum,hvrt) ;
%tricost(vnew,enew,tnew,tnum,hnew) ;

%drawnow;
    
%set(figure(1),'units','normalized', ...
%  'position',[.05,.50,.30,.35]);
%set(figure(2),'units','normalized', ...
%  'position',[.35,.50,.30,.35]);
%set(figure(3),'units','normalized', ...
%  'position',[.05,.05,.30,.35]);
%set(figure(4),'units','normalized', ...
%  'position',[.35,.05,.30,.35]);
