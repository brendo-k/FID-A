function op_CSIPlotVoxelSpec(in, x, y, ppmmin, ppmmax)
    if in.flags.spectralFT == 0 && in.flags.spatialFT == 0
        error('please fourier transfrom along the spatial and spectral dimensions')
    end
    
    if ~exist('ppmmin', 'var')
        ppmmin = min(in.ppm); 
    end
    if ~exist('ppmmax', 'var')
        ppmmax = max(in.ppm);
    end
    
    signal = permute(in.specs, nonzeros([in.dims.t, in.dims.x, in.dims.y, in.dims.coils, in.dims.averages]));
    ppm = in.ppm;
    range = ppm > ppmmin & ppm < ppmmax;
    plot_signal = signal(range, x, y);
    ppm_range = ppm(range);
    plot(ppm_range, real(plot_signal))
end