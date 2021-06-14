function out = op_CSIRemoveWater(in, xlim, ylim)
residual_error = 0;
    for x = 1:in.sz(in.dims.x)
        for y = 1:in.sz(in.dims.y)
            voxel = CSItoMRS(in, x, y);
            voxel_supressed = op_removeWater(voxel, [4.4, 5], 20, 1500, 0);
            in.specs(:, x, y) = voxel_supressed.specs;
            residual_error = residual_error + voxel_supressed.watersupp.residual_error;
        end
    end
    out = in;
    out.watersupp.residual_error = residual_error;

end