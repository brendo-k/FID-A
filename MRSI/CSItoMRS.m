% CSItoTwix.m
% converts a CSI twix to a single voxel FID-A object wich allows for the single 
% voxel FID-A porcessing tools to be used
% 
% USAGE:
% out = CSItoMRS(in, xCoordinate, yCoordinate)
%
% INPUT: 
% in = FID-A CSI object
% xCoordinate = X coordinate wanted to be converted
% yCoordinate = Y coorindate wanted to be converted
%
% OUTPUT:
% out = single voxel FID-A object

function out = CSItoMRS(in, xCoordinate, yCoordinate)
    if(in.flags.spatialFT == 0 || in.flags.spectralFT == 0)
        error('please fourier transform along the spatial and spectral dimension before using this function');
    end
    out = in;
    if(in.dims.coils ~= 0)
        %change the specs to a single voxel with all the coils
        out.specs = in.specs(:,:,xCoordinate, yCoordinate);
    else
        %change the specs to a single voxel
        out.specs = in.specs(:,xCoordinate, yCoordinate);
    end
    out.specs = squeeze(out.specs);
    
    %inverse fourrier transform to get fids at x and y coordiantes
    out.fids = squeeze(ifft(ifftshift(out.specs, in.dims.t), [], in.dims.t));
    %adjust structure to match FID-A object
    out.sz = size(out.fids);
    out.dims.x = 0;
    out.dims.y = 0;
    
end
