classdef Trajectory
    properties
        k_trajectory
        imageSize 
        FoV
        repetitionTime = 1
        pixel_width
        delta_K  
        K_coordinates
        spacial_coordinates
        FovK
        scanTime
        dwellTime
        sw
    end
    methods
        function obj = Trajectory(k_trajectory, imageSize, FoV, dwellTime, ...
            k_dim, sw)
       
            obj.k_trajectory  = k_trajectory;
            obj.pixel_width.x = FoV(1)/imageSize(1);
            obj.pixel_width.y = FoV(2)/imageSize(2);
            obj.imageSize = imageSize;
            obj.FoV.y = FoV(2);
            obj.FoV.x = FoV(1);
            obj.delta_K.x = 1/FoV(1);
            obj.delta_K.y = 1/FoV(2);
            obj.FovK.x = 1/obj.pixel_width.x;
            obj.FovK.y = 1/obj.pixel_width.y;
            obj.dwellTime = dwellTime;
            obj.K_coordinates.k{1} = -obj.FovK.x/2+obj.delta_K.x/2:obj.delta_K.x:obj.FovK.x/2-obj.delta_K.x/2;
            obj.K_coordinates.k{2} = -obj.FovK.y/2+obj.delta_K.y/2:obj.delta_K.y:obj.FovK.y/2-obj.delta_K.y/2;
            if(isnan(obj.K_coordinates.k{1}))
                obj.K_coordinates.k{1} = [0];
            end
            if(isnan(obj.K_coordinates.k{2}))
                obj.K_coordinates.k{2} = [0];
            end
            obj.K_coordinates.dims = k_dim;
            obj.spacial_coordinates.x = -obj.FoV.x/2+obj.pixel_width.x/2:obj.pixel_width.x:obj.FoV.x/2-obj.pixel_width.x/2;
            obj.spacial_coordinates.y = -obj.FoV.y/2+obj.pixel_width.y/2:obj.pixel_width.y:obj.FoV.y/2-obj.pixel_width.y/2;
            obj.scanTime = size(k_trajectory,1)*(obj.repetitionTime + size(k_trajectory,2)*dwellTime);
            obj.sw = sw;
        end
    end
end