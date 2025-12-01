% ========================================================================
% 3D LUNG MODEL + TUMOR SNAPSHOTS + GIF (CLEAN FIXED VERSION)
% ========================================================================
clear; close all; clc;

% ===================== CREATE 3D LUNG GEOMETRY ========================

[x,y,z] = meshgrid( ...
    linspace(-11, 11, 100), ...
    linspace(-7, 7, 100), ...
    linspace(-24, 0, 100));

F_right = ((x-2)/3.5).^2 + (y/2.5).^2 + ((z+7)/6).^2 - 1;
mask_right = F_right <= 0;

F_left = ((x+2)/3.2).^2 + (y/2.3).^2 + ((z+7)/6).^2 - 1;
mask_left = F_left <= 0;

mask_lungs = mask_right | mask_left;
lung_voxels = find(mask_lungs);
Nvox = length(lung_voxels);

lungs_val = double(mask_lungs);

% ========================== LOAD CSV DATA =============================
healthy = readmatrix('healthy_results.csv');
smoker  = readmatrix('smoker_results.csv');

t_h = healthy(:,1); 
c_h = healthy(:,2);

t_s = smoker(:,1);
c_s = smoker(:,2);

% ========================= FOLDER + GIF SETUP ==========================
output_folder = "frames_lungs";
if ~exist(output_folder, "dir")
    mkdir(output_folder);
end

gif_name = "tumor_evolution.gif";
if exist(gif_name,'file')
    delete(gif_name);
end

% ======================= SCALING FACTOR ================================
scale = 14;   % 1 voxel = 100 cancer cells

% ======================== GENERATE FRAMES ==============================
fprintf("Generating images and GIF...\n");
counter = 1;

% --------------------------- HEALTHY ------------------------------------
for step = 1:length(t_h)

    Ntumor = round(c_h(step) / scale);
    Ntumor = min(Nvox, max(0, Ntumor));

    save_frame(x,y,z,lungs_val,lung_voxels,Ntumor,step,t_h(step),'HEALTHY',output_folder);

    img = imread(sprintf('%s/frame_HEALTHY_%04d.png', output_folder, step));
    [A,map] = rgb2ind(img,256);

    if counter == 1
        imwrite(A,map,gif_name,'gif','LoopCount',inf,'DelayTime',0.18);
    else
        imwrite(A,map,gif_name,'gif','WriteMode','append','DelayTime',0.18);
    end

    counter = counter + 1;
end

% --------------------------- SMOKER ------------------------------------
for step = 1:length(t_s)

    Ntumor = round(c_s(step) / scale);
    Ntumor = min(Nvox, max(0, Ntumor));

    save_frame(x,y,z,lungs_val,lung_voxels,Ntumor,step,t_s(step),'SMOKER',output_folder);

    img = imread(sprintf('%s/frame_SMOKER_%04d.png', output_folder, step));
    [A,map] = rgb2ind(img,256);
    imwrite(A,map,gif_name,'gif','WriteMode','append','DelayTime',0.18);
end

% ======================== SAVE FRAME FUNCTION ==========================
function save_frame(x,y,z,lungs_val,lung_voxels,Ntumor,step,t,who,folder)

    Tumor = zeros(size(lungs_val));

    if Ntumor > 0
        idx = lung_voxels(1:Ntumor);   % deterministic
        Tumor(idx) = 1;
    end

    fv_tumor = isosurface(x,y,z,Tumor,0.5);

    fig = figure('Visible','off');
    hold on;

    fv_l = isosurface(x,y,z,lungs_val,0.5);
    patch(fv_l,'FaceColor',[1 0.7 0.7],'EdgeColor','none','FaceAlpha',0.2);

    if ~isempty(fv_tumor.vertices)
        patch(fv_tumor,'FaceColor',[1 0 0],'EdgeColor','none','FaceAlpha',0.9);
    end

    daspect([1 1 1]);
    view(40,20);
    camlight headlight;

    title(sprintf('%s – t = %.1f days — %d tumor voxels', who, t, Ntumor));
    axis off;

    filename = sprintf('%s/frame_%s_%04d.png', folder, who, step);
    saveas(fig, filename);
    close(fig);
end
