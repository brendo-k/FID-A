%Simulating the trjaectory of an EPSI trajectory. Only simulating EPSI
%trajectory that aquires in one direction (no aquasition during rewind
%gradients).
%Input: par (optional)
%           consists of fields:
%               dwellTime: dwell time in seconds
%               Fov: desired fov size in meters
%               imageSize: desired resolution. Array with 3 entries corresponding to x, y,and spectral axises respectivly
%               readOutTime: readout time in the spectral dimension
%               repetitionTime: Time between each excitation
%
%Output:    Traj: 2d matrix of the resulting K-space trajectory from the input
%              params. Dimensions are excitation number and readout for the frist and
%              second dimensions respectivly.
%           scanTime: Resulting scan time from trajectory. In seconds
%           par: paramaters (same fields as above) used for EPSI simulation   

function [traj, scanTime, par] = EPSI(par)

    if(nargin < 1)
        par.dwellTime = 31.25e-6; %[s]
        par.Fov = [0.18, 0.18]; %FoV in m in the x and y directions
        par.imageSize = [32, 32, 40]; %voxels resolution in the x and y and spectral dimension
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
    
    %initalize empty array for trajectory 
    traj = zeros(par.imageSize(2), par.imageSize(1)*par.spectralPoints);
    
    %readout along the first dimension
    for j = 1:par.imageSize(2)
        %multiple reads of the x axis
        traj(j,:) = repmat(kSpaceX, [1,par.imageSize(3)]);
        
        %add the corresponding y position
        traj(j,:) = traj(j,:) + 1i*kSpaceY(j);
    end
    
    %Calculate scan time.
    %(dwelltime * imagesize * 2) is the time to get back to the same position
    %in k space which needs to be repeated equivalently to spectral points.
    scanTime = (dwellTime*par.imageSize(1)*2*par.spectralPoints+par.repetitionTime)*par.imageSize(2);

end
    
