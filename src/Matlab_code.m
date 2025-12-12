clear; close all; clc;

% First, let's create the 3d lung geometry 

[x,y,z] = meshgrid( ...
    linspace(-11, 11, 100), ... % left-right length: 22cm corresponds to an average man lung 
    linspace(-7, 7, 100), ... % front-back length: 14cm  
    linspace(-24, 0, 100)); % % up-down length: 24cm 
% the term "100" corresponds to the number of points along each line 
    

F_right = ((x-2)/3.5).^2 + (y/2.5).^2 + ((z+7)/6).^2 - 1; % making the ellipsoid for each lung 
mask_right = F_right <= 0; % defines which points of the domain are considered "the lung"; <=0 is inside, giving the cells inside the attribute "true" and the ones outside the attribute "false" 

F_left = ((x+2)/3.2).^2 + (y/2.3).^2 + ((z+7)/6).^2 - 1; 
mask_left = F_left <= 0;

mask_lungs = mask_right | mask_left; % combine the two ellipsoids 
lung_voxels = find(mask_lungs); % gets all indices of the points that are in the mask of the lungs 
Nvox = length(lung_voxels); % counts the number of lung voxels 

lungs_val = double(mask_lungs); % converts the "true" and "false" in the mask to 1 and 0, numerical values are needed for making the isosurface next 

% loading the csv data obtained with the C code 

healthy = readmatrix('healthy_results.csv');
smoker  = readmatrix('smoker_results.csv');

t_h = healthy(:,1); % each time step for the healthy lung
c_h = healthy(:,2); % each cell count for the healthy lung 

t_s = smoker(:,1);
c_s = smoker(:,2);

% output folder and GIF making 

output_folder = "frames_lungs";
if ~exist(output_folder, "dir")
    mkdir(output_folder);
end 
% if the folder does not exist, it creates it and if it does, it doesn't do anything 


gif_name = "tumor_evolution.gif";
if exist(gif_name,'file')
    delete(gif_name);
end 
% if there already is a gif in the folder it deletes it and makes a new one

% scaling factor for visualisation, so the lung doesn't get full too soon but still has a cancer progress that can be seen clearly 

scale = 14;   % 1 voxel = 100 cancer cells

% generating the frames 

counter = 1; % variable to keep track of which frame we are on, starting with the first 

% healthy case 

for step = 1:length(t_h) % for each time step, or line in our csv file 

    Ntumor = round(c_h(step) / scale); % counts the number of tumorous voxels ( number of cancerous cells at the time step divided by the scale ) 
    Ntumor = min(Nvox, max(0, Ntumor)); % keeps our Ntumor in a reasonable range, not more than the actual number of voxels in the lung but also not a negative number of cancerous cells 

    save_frame(x,y,z,lungs_val,lung_voxels,Ntumor,step,t_h(step),'healthy',output_folder); % makes an 'image' with the state of the cancer at this time step 

    img = imread(sprintf('%s/frame_healthy_%04d.png', output_folder, step)); % reads the PNG image we just created 
    [A,map] = rgb2ind(img,256); % converts the RGB image to index colors 

    if counter == 1
        imwrite(A,map,gif_name,'gif','LoopCount',inf,'DelayTime',0.24); % creates the GIF, makes it loop forever (inf), pause 0.18 seconds between frames, the first one 'builds' the gif
    else
        imwrite(A,map,gif_name,'gif','WriteMode','append','DelayTime',0.24); % appends each image to the first one to create the whole GIF 
    end

    counter = counter + 1; % indicates we are moving to the next frame 
end

% smoker case, same structure of code as the one for healthy

for step = 1:length(t_s)

    Ntumor = round(c_s(step) / scale);
    Ntumor = min(Nvox, max(0, Ntumor));

    save_frame(x,y,z,lungs_val,lung_voxels,Ntumor,step,t_s(step),'smoker',output_folder);

    img = imread(sprintf('%s/frame_smoker_%04d.png', output_folder, step));
    [A,map] = rgb2ind(img,256);
    imwrite(A,map,gif_name,'gif','WriteMode','append','DelayTime',0.30);
end

% save frame function 

function save_frame(x,y,z,lungs_val,lung_voxels,Ntumor,step,t,who,folder) % who = string, 'healthy' or 'smoker' 

    Tumor = zeros(size(lungs_val)); % creates an array of 0  and 1, this will become the tumor mask, 0 = no tumor, 1 = tumor voxel  

    if Ntumor > 0 % if the tumor has at least one voxel 
        idx = lung_voxels(1:Ntumor);   % take the first Ntumor lung voxel positions
        Tumor(idx) = 1; % set them to one in the tumor mask, so mark them as tumorous 
    end

    fv_tumor = isosurface(x,y,z,Tumor,0.5); % creates a 3D surface mesh around the tumor volume

    fig = figure('Visible','off'); % creates a figure without showing it
    hold on; % allows to add to add more surfaces on the figure 

    fv_l = isosurface(x,y,z,lungs_val,0.5); % 3d mesh of the lung 
    patch(fv_l,'FaceColor',[1 0.7 0.7],'EdgeColor','none','FaceAlpha',0.2); % draws it light pink with a transparency of 20 percent  

    if ~isempty(fv_tumor.vertices)
        patch(fv_tumor,'FaceColor',[1 0 0],'EdgeColor','none','FaceAlpha',0.9);
    end 
    % if the tumor mesh is not empty, draw it in red 

    daspect([1 1 1]);
    view(40,20);
    camlight headlight;
    % 3d appearance 
    
    title(sprintf('%s – t = %.1f days — %d tumor voxels', who, t, Ntumor));
    axis off; % removes axis for a better visualisation 

    filename = sprintf('%s/frame_%s_%04d.png', folder, who, step); % filename for saving 
    saveas(fig, filename); % writes the png file 
    close(fig); % deletes the figure from memory 
end




