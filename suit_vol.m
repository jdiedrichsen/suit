function V = suit_vol(source,varargin)
% Calculate ROI volume in one or more images in native space, returns  
% the number of voxels as well as the mm3 based in the information in the 
% affine matrix. Always check that the calculated voxel volume is correct.
%
% Input:    Cell array with binary images(Masks) to calculate volume or
%           resliced cerebellar atlas. No need of cell for single image.
%
% Atlas:    If you use the flag 'Atlas' as a second argument, the function
%           will return the volume of each of the regions in the
%           SUIT atlas. The Atlas must be resliced into the native space
%           using suit_normalize_dartel_inv. For stats calculation 
%           in SUIT space use suit_lobuli_summarize.
%
% Output:   Structure including the following fields:
%       name:   Filename of the source
%       dir:    Directory of the source
%       vox:    Number of voxels per region
%       vmm:    Volume in mm3 per region
%       Vsize:  Volume of one voxel
%       reg:    List of regions of SUIT atlas
%
% Example:
%
%   Cerebellar volume using the mask from suit_isolate_seg:
%
%       V = suit_vol({'c_<name>_pcereb.nii'});   
%
%   Volume of each lobuli using the resliced atlas in native space:
%   
%       V = suit_vol({'iw_Cerebellar-SUIT-<name>'},'Atlas');
% _______________________________________________________________________
% Carlos R. Hernandez-Castillo 2018

if (iscell(source)); source=char(source); end;

if isempty(varargin)
    for i = 1:size(source,1)
        input = source(i,:);
        [source_dir,Sname,~,~]=spm_fileparts(input);
        if (isempty(source_dir))
            source_dir=pwd;
        end;

        V1 = spm_vol(input);
        VoxSize = abs(V1.mat(1,1) * V1.mat(2,2) * V1.mat(3,3));
        v1 = spm_read_vols(V1);
        vox = numel(find(single(v1)==1));
        vmm = vox * VoxSize;

        V(i).name = Sname;
        V(i).dir  = source_dir;
        V(i).vox  = vox;
        V(i).vmm  = vmm;
        V(i).Vsize= VoxSize;
    end
    return
    
end

namereg={'Left I_IV','Right I_IV','Left V','Right V',...
    'Left VI','Vermis VI','Right VI',...
    'Left Crus I','Vermis Crus I','Right Crus I',...
    'Left Crus II','Vermis Crus II','Right Crus II',...
    'Left VIIb','Vermis VIIb','Right VIIb',...
    'Left VIIIa','Vermis VIIIa','Right VIIIa',...
    'Left VIIIb','Vermis VIIIb','Right VIIIb',...
    'Left IX','Vermis IX','Right IX',...
    'Left X','Vermis X','Right X','Left Dentate',...
    'Right Dentate','Left Interposed','Right Interposed',...
    'Left Fastigial','Right Fastigial'};

if strcmp(varargin{1},'Atlas')
    for i = 1:size(source,1)
        input = source(i,:);
        [source_dir,Sname,~,~]=spm_fileparts(input);
        if (isempty(source_dir))
            source_dir=pwd;
        end;
        V1 = spm_vol(input);
        VoxSize = abs(V1.mat(1,1) * V1.mat(2,2) * V1.mat(3,3));
        v1 = spm_read_vols(V1);

        for r=1:34
            nvox(r) = length(find(v1==r));
            vmm(r) = nvox(r) * VoxSize;
        end
        
        V(i).name = Sname;
        V(i).dir  = source_dir;
        V(i).vox  = nvox;
        V(i).vmm  = vmm;
        V(i).Vsize= VoxSize;
        V(i).reg  = namereg;

    end
    
end


