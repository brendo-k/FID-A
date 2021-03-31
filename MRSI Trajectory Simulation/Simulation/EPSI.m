%Simulating the trjaectory of an EPSI trajectory. Only simulating EPSI
%trajectory that aquires in one direction (no aquasition during rewind
%gradients).
%Input: par (optional)
%           consists of fields:
%             dwellTime: dwell time in seconds
%             Fov: desired fov size in meters
%             imageSize: desired resolution. Array with 2 entries corresponding to x and y axises respectivly
%             readOutTime: readout time in the spectral dimension
%             repetitionTime: Time between each excitation

%TODO: Calculate ramping time to each k space position to add to total scan
%time
            

function [traj, scanTime, par] = EPSI(par)

    if(nargin < 1)
        par.dwellTime = 31.25e-6; %[s]
        par.Fov = [0.18, 0.18]; %FoV in m in the x and y directions
        par.imageSize = [32, 32]; %voxels resolution in the x and y direction
        par.spectralPoints = 40; % Number of readout points in the spectral dimension
        par.repetitionTime = 1; %s
    end
    %calculating the same image parameters as the default parameters in Rosette.m 

    %updating local variables from parameters
    Fov = par.Fov;
    deltaFovX = Fov(1)/par.imageSize(1); %[m]
    deltaFovY = Fov(2)/par.imageSize(2); %[m]
    deltaKX = 1/Fov(1); %[1/m]
    deltaKY = 1/Fov(2); %[1/m]
    FovKX = 1/deltaFovX; %[1/m]
    FovKY = 1/deltaFovY; %[1/m]
    dwellTime = par.dwellTime; %[s]
    sw = 1/(dwellTime*par.imageSize(1)*2); %[Hz]
    
    
    %calculating trajectory for each shot
    kSpaceX = -FovKX/2 + deltaKX/2:deltaKX:FovKX/2 - deltaKX/2;
    kSpaceY = -FovKY/2 + deltaKY/2:deltaKY:FovKY/2 - deltaKY/2;
    
    traj = zeros(par.imageSize(2), par.imageSize(1)*par.spectralPoints);
    
    for j = 1:par.imageSize(2)
        traj(j,:) = repmat(kSpaceX, [1,par.spectralPoints]) + 1i*kSpaceY(j);
    end
    
    scanTime = (dwellTime*par.imageSize(1)*2*par.spectralPoints+par.repetitionTime)*par.imageSize(2);

end
    
