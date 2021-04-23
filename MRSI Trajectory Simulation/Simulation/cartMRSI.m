%Simulating the trjaectory of the basic cartesian MRSI trajectory
%Input: par (optional)
%           consists of fields:
%             dwellTime: dwell time in seconds
%             Fov: desired fov size in meters
%             imageSize: desired resolution. Array with 2 entries corresponding to x and y axises respectivly
%             readOutTime: readout time in the spectral dimension
%             repetitionTime: Time between each excitation

%TODO: Calculate ramping time to each k space position to add to total scan
%time
            

function [obj] = cartMRSI(par)

    if(nargin < 1)
        par.dwellTime = 125/100000; %us ->s
        par.Fov = [0.18, 0.18]; %FoV in m
        par.imageSize = [7, 7, 40]; %voxels in the x and y direction
    end
    %calculating the same image parameters as the default parameters in Rosette.m 

    %updating local variables from parameters
    Fov = par.Fov;
    deltaFovX = Fov(1)/par.imageSize(1);
    deltaFovY = Fov(2)/par.imageSize(2);
    deltaKX = 1/Fov(1);
    deltaKY = 1/Fov(2);
    FovKX = 1/deltaFovX;
    FovKY = 1/deltaFovY;
    readOutTime = par.dwellTime*(par.imageSize(3)-1);
    dwellTime = par.dwellTime;
    
    %calculating trajectory for each shot
    kSpaceX = -FovKX/2+deltaKX/2:deltaKX:FovKX/2-deltaKX/2;
    kSpaceY = -FovKY/2+deltaKY/2:deltaKY:FovKY/2-deltaKY/2;
    [meshx, meshy] = meshgrid(kSpaceX, kSpaceY);
    meshy = meshy .* 1i;
    traj = meshx + meshy;
    
    dim = size(traj);
    excite = dim(1) * dim(2);
    
    traj = reshape(traj, [excite ,1]);
    
    %add plus one because the zeroth point is counted
    traj = repmat(traj, [1,readOutTime/dwellTime + 1]);
    
    obj = Trajectory(traj, par.imageSize, Fov, readOutTime, ["x", "y"], 1/dwellTime);
end
    
