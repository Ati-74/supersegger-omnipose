function data = intMakeRegs( maskpath, data, CONST, mask_bad_regs, good_regs )
% intMakeRegs : creates info for bad regions or makes new regions
%
% INPUT :
%       data : cell file (seg/err file)
%       CONST : segmentation constants
%       mask_bad_regs : mask of bad regions (their score is set to 0)
%       good_regs : if 1 all scores are set to 1
% OUTPUT : 
%       data : cell file with region fields
%
% Copyright (C) 2016 Wiggins Lab 
% Written by Stella Stylianidou & Paul Wiggins.
% University of Washington, 2016
% This file is part of SuperSegger.
% 
% SuperSegger is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% SuperSegger is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with SuperSegger.  If not, see <http://www.gnu.org/licenses/>.

E = CONST.regionScoreFun.E;

%sets all scores to 1
if ~exist('good_regs','var') || isempty(good_regs)
    good_regs = false;
end

%ss = size( data.mask_cell ); %comment out
NUM_INFO = CONST.regionScoreFun.NUM_INFO;
%%

% Construct the Raw Image Path
[maskDir, maskName, ~] = fileparts(maskpath); % Extract directory and name of mask

% List all the PNG files in the directory
all_files = dir(fullfile(maskDir, '*.png'));
all_filenames = {all_files.name};

% Find the index of the specified image in the list
file_index = find(strcmp(all_filenames, strcat(maskName, '.png')));

[parentDir, ~, ~] = fileparts(maskDir);  % Get the parent directory
[parentDir, ~, ~] = fileparts(parentDir);  % Get the parent directory
rawImageDir = fullfile(parentDir, 'raw_im'); % Get directory of raw images

% crop box
load(strcat(rawImageDir, '/cropbox.mat'),'crop_box_array');

labeledmask = double(imread(maskpath)); %load cellpose mask

stats = regionprops(labeledmask, 'PixelList'); % Get pixel coordinates of each region

for k = 1:length(stats)
    x = stats(k).PixelList(:,1);
    y = stats(k).PixelList(:,2);
    
    if any(x <= round(crop_box_array{1, 1}(file_index , 2))) || ...
            any(x >= round(crop_box_array{1, 1}(file_index , 4)) - 1) || ...
            any(y <= round(crop_box_array{1, 1}(file_index , 1))) || ...
            any(y >= round(crop_box_array{1, 1}(file_index , 3)) - 1)
   
        % Populate the mask with the region
        diff_x_values = diff(x);
        if any(diff_x_values > 2)
            if any(x(find(diff_x_values > 2) + 1:length(x)) <= round(crop_box_array{1, 1}(file_index , 2)) + 1) || ...
                    any(x(find(diff_x_values > 2) + 1:length(x)) >= round(crop_box_array{1, 1}(file_index , 4)) - 1) || ...
                    any(y(find(diff_x_values > 2) + 1:length(x)) <= round(crop_box_array{1, 1}(file_index , 1)) + 1) || ...
                    any(y(find(diff_x_values > 2) + 1:length(x)) >= round(crop_box_array{1, 1}(file_index , 3)) - 1)
                
                for j = find(diff_x_values > 2) + 1: length(x)
                    labeledmask(y(j), x(j)) = 0;
                end  
            end
            
            if any(x(1:find(diff_x_values > 2)) <= round(crop_box_array{1, 1}(file_index , 2)) + 1) || ...
                    any(x(1:find(diff_x_values > 2)) >= round(crop_box_array{1, 1}(file_index , 4)) - 1) || ...
                    any(y(1:find(diff_x_values > 2)) <= round(crop_box_array{1, 1}(file_index , 1)) + 1) || ...
                    any(y(1:find(diff_x_values > 2)) >= round(crop_box_array{1, 1}(file_index , 3)) - 1)
                
                for j = 1:find(diff_x_values > 2)
                    labeledmask(y(j), x(j)) = 0;
                end  
            end
        else
           for j = 1:length(x)
                labeledmask(y(j), x(j)) = 0;
           end 
        end
    end
end


reduced_labeledmask = labeledmask(round(crop_box_array{1, 1}(file_index , 1)) : round(crop_box_array{1, 1}(file_index , 3)) - 1,  round(crop_box_array{1, 1}(file_index , 2)): round(crop_box_array{1, 1}(file_index , 4)) - 1);

%% my script to write the segmented objects
% Generate unique colors for each region
unique_colors = generate_unique_colors(max(reduced_labeledmask(:)));

% Create an RGB image based on the unique colors
[rows, cols] = size(reduced_labeledmask);
colored_mask = zeros(rows, cols, 3); % Initialize an RGB image

for k = 1:max(reduced_labeledmask(:))
    mask = reduced_labeledmask == k;
    for channel = 1:3
        colored_mask(:, :, channel) = colored_mask(:, :, channel) + mask * unique_colors(k, channel);
    end
end
colored_mask = uint8(colored_mask); % Convert to ui


% Determine directory and filename for saving
[out_dir, filename, ~] = fileparts(maskpath);
[parent_dir, ~, ~] = fileparts(out_dir);

save_dir = fullfile(parent_dir, 'rgbObjects');
save_dir_BW = fullfile(parent_dir, 'BWObjects');


% Create directory save_dir if it doesn't exist
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end


if ~exist(save_dir_BW, 'dir')
    mkdir(save_dir_BW);
end


% Save the result in .mat format
save_path = fullfile(save_dir, [filename '.mat']);
save(save_path, 'colored_mask');

%% Create Binary BW Image
bw_image = reduced_labeledmask > 0; % Any non-zero value in the mask is considered as object

% Save the binary image as .tif
bw_save_path = fullfile(save_dir_BW, [filename '.tif']);
imwrite(bw_image, bw_save_path);

ss = size(labeledmask);
data.regs.regs_label = labeledmask; %create regs_label with cellpose mask 

%data.regs.regs_label = bwlabel( data.mask_cell );
data.regs.num_regs = max( data.regs.regs_label(:) );
data.regs.props = regionprops( data.regs.regs_label, ...
    'BoundingBox','Orientation','Centroid','Area');
data.regs.score  = ones( data.regs.num_regs, 1 );
data.regs.scoreRaw = ones( data.regs.num_regs, 1 );
data.regs.info = zeros( data.regs.num_regs, NUM_INFO );


for ii = 1:data.regs.num_regs
    
    [xx,yy] = getBBpad( data.regs.props(ii).BoundingBox, ss, 1);
    mask = data.regs.regs_label(yy,xx)==ii;
    data.regs.info(ii,:) = CONST.regionScoreFun.props( mask, data.regs.props(ii) );
    
    if exist( 'mask_bad_regs', 'var' ) && ~isempty( mask_bad_regs )
        data.regs.scoreRaw(ii) = CONST.regionScoreFun.fun(data.regs.info(ii,:), E);
        data.regs.score(ii) = data.regs.scoreRaw(ii) > 0;
        mask_ = mask_bad_regs(yy,xx);
        if any( mask(mask_) )
            data.regs.score(ii) = 0;
        end
    end
    
end

if ~exist( 'mask_bad_regs', 'var' ) || isempty( mask_bad_regs )
    if good_regs
        data.regs.scoreRaw = CONST.regionScoreFun.fun(data.regs.info, E)';
        data.regs.score = ones( data.regs.num_regs, 1 );
    else
        data.regs.scoreRaw = CONST.regionScoreFun.fun(data.regs.info, E)';
        data.regs.score =   data.regs.scoreRaw>0;
    end
end



% Function to generate unique colors
function colors = generate_unique_colors(num_colors, max_rgb)
    if nargin < 2
        max_rgb = 240;
    end

    skip = ceil(max_rgb / nthroot(num_colors, 3));

    colors = zeros(num_colors, 3); % Initialize an Nx3 matrix for RGB values
    color_count = 1;

    while color_count <= num_colors
        for r = 0:skip:max_rgb
            for g = 0:skip:max_rgb
                for b = 0:skip:max_rgb
                    if color_count <= num_colors
                        colors(color_count, :) = [r+1, g+1, b+1];
                        color_count = color_count + 1;
                    else
                        return;
                    end
                end
            end
        end

        skip = skip - 1;
        if skip <= 0
            error('Unable to generate the required number of colors with the given constraints.');
        end
    end
end

end
