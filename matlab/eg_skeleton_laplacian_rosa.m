function P = eg_skeleton_laplacian_rosa(filename)
% extract curve skeleton from a point cloud or triangular mesh
% update: 2010-8-19
% update: 2010-7-12
% create: 2009-4-26
% by: JJCAO, deepfish @ DUT
%
%% setting
addpath(genpath(fullfile(fileparts([mfilename('fullpath'), '.m']), 'toolbox')))
options.USING_POINT_RING = GS.USING_POINT_RING;
extension='.off';

%% Step 0: read file (point cloud & local feature size if possible), and
% normalize the modle.

tic
P.filename = [filename extension];% point set
[P.pts,P.faces] = read_mesh(P.filename);
P.npts = size(P.pts,1);
if exist([filename '_fe.txt'],'file') % result of Tamal K Dey's NormFet
    P.radis = load([filename '_fe.txt']);
else
    P.radis = ones(P.npts,1);
end

P.pts = GS.normalize(P.pts);
[P.bbox, P.diameter] = GS.compute_bbox(P.pts);
fprintf('read point set:\n');
toc

%% Step 1: build local 1-ring
% build neighborhood, knn?
tic
P.k_knn = GS.compute_k_knn(P.npts);
if options.USING_POINT_RING
    P.rings = compute_point_point_ring(P.pts, P.k_knn, []);
else    
    P.frings = compute_vertex_face_ring(P.faces);
    P.rings = compute_vertex_ring(P.faces, P.frings);
end
fprintf('compute local 1-ring:\n');
toc

%% Step 1: Contract point cloud by Laplacian
tic
[P.cpts, t, initWL, WC, sl] = contraction_by_mesh_laplacian(P, options); %#ok<ASGLU>
fprintf('Contraction:\n');
toc
%% Step 2: Point to curve C by cluster ROSA2.0
tic
P.sample_radius = P.diameter*0.02;
P = rosa_lineextract(P,P.sample_radius, 1);
fprintf('to curve:\n');
toc
%% Show results
figName='Original point cloud and its contraction';
figure('Name',figName,'Color','w','Renderer','OpenGL'); movegui('onscreen')
scatter3(P.pts(:,1),P.pts(:,2),P.pts(:,3),30,'.','MarkerEdgeColor', GS.PC_COLOR);  hold on;
scatter3(P.cpts(:,1),P.cpts(:,2),P.cpts(:,3),30,'.r');
axis off;axis equal;

figName='Original point cloud and its skeleton';
figure('Name',figName,'Color','w','Renderer','OpenGL'); movegui('onscreen')
scatter3(P.pts(:,1),P.pts(:,2),P.pts(:,3),20,'.','MarkerEdgeColor', GS.PC_COLOR);  hold on;
showoptions.sizep=400;showoptions.sizee=2;
plot_skeleton(P.spls, P.spls_adj, showoptions);
axis off;axis equal;