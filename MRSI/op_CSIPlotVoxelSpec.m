function fig = op_CSIPlotVoxelSpec(in, x, y, ppmmin, ppmmax, xlab, ylab, title)
    if in.flags.spectralFT == 0 && in.flags.spatialFT == 0
        error('please fourier transfrom along the spatial and spectral dimensions')
    end
    if ~exist('x', 'var') || ~exist('y', 'var')
        error('please provide an x and y coordinate');
    end
    
    single_vox = CSItoMRS(in, x, y);
    if(~exist('ppmmin', 'var') || ~exist('ppmmax','var'))
        fig = op_plotspec(single_vox);
    elseif(~exist('xlab', 'var'))
        fig = op_plotspec(single_vox, ppmmin, ppmmax);
    elseif(~exist('ylab', 'var'))
        fig = op_plotspec(single_vox, ppmmin, ppmmax, xlab);
    elseif(~exist('title', 'var'))
        fig = op_plotspec(single_vox, ppmmin, ppmmax, xlab, ylab);
    else
        fig = op_plotspec(single_vox, ppmmin, ppmmax, xlab, ylab, title);
    end
    
end