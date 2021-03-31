function [TrajectoryStruct] = CreateTrajectory(Trajectory, par)
    
    if(lower(Trajectory) == "cartesian")
        %Using cartMRSI to generate sructure
        [traj, scanTime] = cartMRSI(par);
        TrajectoryStruct.scan_time = scanTime;
        TrajectoryStruct.trajectory = traj; 
    
    elseif (lower(Trajectory) == "rosette")
        %Using Rosette fucntion to generate strucutre
        [~, gradient, finalKSpaceTraj, gMax, scanTime, maxSlewRate, rampingPoints] = Rosette(par);
        TrajectoryStruct.gradient_trajectory = gradient(:);
        TrajectoryStruct.k_space_trajectory = finalKSpaceTraj(:);
        TrajectoryStruct.scan_time = scanTime;
        TrajectoryStruct.max_gradient = gMax;
        TrajectoryStruct.max_slew_rate = maxSlewRate;
        
        %Field of view
        TrajectoryStruct.FoV = 1/sqrt((finalKSpaceTraj(1,1) - finalKSpaceTraj(1,2))^2 ...
                                + (finalKSpaceTraj(2,1) - finalKSpaceTraj(2,2))^2);
        %voxel size
        TrajectoryStruct.voxel_size = 1/par.kMax;
        TrajectoryStruct.total_points = size(finalKSpaceTraj,1);/?fbclid=IwAR0q-FAu3O-nlsBTvStKzXvacMGRTqu9b2SYOIL4brg2_xb0NeQ7A_6IDus
        
        %number of ramping points
        TrajectoryStruct.ramping_points = rampingPoints;
        
        %number of shots
        TrajectoryStruct.number_of_shots = finalKSpaceTraj.size(2);
        
        %time for one segment
        %time for the trajectory to return back to the original position
        %found by the equation K_max*sin(omega1*t)e^(i*oega2*t)
        %time between spectral sampling
        TrajectoryStruct.delat_t = pi/par.omega1;
        
        %number of points per segment
        TrajectoryStruct.points_per_segment = size(0:par.dwellTime:TrajectoryStruct.time_for_segment, 1);
        
        %number of points for spectral readout
        TrajectoryStruct.spectral_readout_poits = size(finalKSpaceTraj,1)/TrajectoryStruct.points_per_segment;
        
        %Dwell time 
        %time between spatial sampling
        TrajectoryStruct.dwell_time = par.dwellTime;
        
        %number of points per fid
        TrajectoryStruct.points_per_fid = size(finalKSpaceTraj, 1);
    
    elseif(~isstring(Trajectory) && ismatrix(Trajectory))
        TrajectoryStruct.k_space_trajectory = Trajectory;
        gyromagneticRatio = 2.675222005e8;
        gradient_x = diff(Trajectory(:,1)/(par.dwell_time * gyromagneticRatio));
        gradient_y = diff(Trajectory(:,2)/(par.dwell_time * gyromagneticRatio));
        TrajectoryStruct.gradient_trajectory(:,1) = gradient_x;
        TrajectoryStruct.gradient_trajectory(:,2) = gradient_y;
        TrajectoryStruct.scan_time = par.scan_time;
    
    else
        error("Please enter in an array for the trajectory or enter the name of a premade trajectory")
    end
        