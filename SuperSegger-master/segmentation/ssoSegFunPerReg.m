function [data, err_flag] = ssoSegFunPerReg( phase, CONST, header, dataname, crop_box)
% ssoSegFunPerReg : starts segmentation of phase image and sets error flags
% It creates the first set of good, bad and permanent segments and if
% CONST.seg.OPTI_FLAG is set to true it optimizes the region sizes.
% It uses perRegionOpti to optimize the regions.
% 
% INPUT :
%       phase_ : phase image
%       CONST : segmentation constants
%       header : string displayed with infromation
%       dataname : 
%       crop_box : information about alignement of the image
% 
%  OUTPUT :
%       data : contains information about the segments and mask, for more
%       information look at superSeggerOpti.
%       err_flag : set to true if there are more segments than max
%       
%
% Copyright (C) 2016 Wiggins Lab 
% Written by Stella Styliandou & Paul Wiggins.
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

if ~exist('header','var')
    header = '';
end

if ~exist('dataname','var')
    dataname = '';
end

if ~exist('crop_box','var')
    crop_box = '';
end

%% skip segmentation/superSeggerOpti
%create the masks and segments
% disp_seg = ~CONST.seg.OPTI_FLAG;
% data = superSeggerOpti( phase ,[], disp_seg ,CONST, 1, header, crop_box);
% 
% if numel(data.segs.score) > CONST.superSeggerOpti.MAX_SEG_NUM;
%     err_flag = true;
%     save([dataname,'_too_many_segs'],'-STRUCT','data');
%     disp( [header,'BSSO ',dataname,'_too_many_segs'] );
%     return
% else
%     err_flag = false;    
% end

%%
% optimize the regions with bad scores
% if CONST.seg.OPTI_FLAG %might need to comment this out since it uses seg?
%     data = perRegionOpti( data, 1, CONST,header);
%     drawnow
% else
    data = [];
    data.phase = phase;
    % masks now stored adjacent to 'phase' directory 
    if exist([extractBefore(dataname,'seg') 'cp_masks' filesep],'dir')
        maskdir = [extractBefore(dataname,'seg') 'cp_masks' filesep];
    elseif exist([extractBefore(dataname,'seg') 'masks' filesep],'dir')
        maskdir = [extractBefore(dataname,'seg') 'masks' filesep];
    end
%     maskdir = [extractBefore(dataname,'seg') 'phase' filesep];
    filename = char(extractBetween(dataname,['seg' filesep],'_seg.mat')); %generalized to filesep 
    maskpath = strcat(maskdir,filename,'c1_cp_masks.png'); %get path of cellpose mask  
    data = intMakeRegs( maskpath, data, CONST ); %input cellpose mask 
    err_flag = false; %shh no errors here...
% end

end

