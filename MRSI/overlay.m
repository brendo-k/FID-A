%overlay.m
% Used as a callback function in sim_CSIoverlayMRI to overaly CSI data onto 
% MRI. The overlay is done with spm_image first, so spm needs to be installed.
% Coordinates are in RAS or neurological coordinates (ie, right = x, anterior = y, superior = z.
%
% USAGE: 
% overlay((optional) plot_type, (optional) in, (optional) coilNum)  
% (must call spm_image before using overlay)
%
% INPUT: 
% plot_type         = a char array of 'real', 'magnitude', and 'imaginary
%                   to be plotted
% in                = MRSI object
% coilNum           = coilNum to plot if coils have not been combined.
%                   (default value)  1 

function overlay(plot_type, in, ppmmin, ppmmax, coilNum)

%call global variable from spm
global st;
%default value of 1 for coilNum

global type;
global csi_obj;
global ppm_min;
global ppm_max;
global coil;
if(exist('plot_type', 'var'))
    type = plot_type;
end
if(exist('in','var'))
    csi_obj = in;
end
if(exist('ppmmin','var'))
    ppm_min = ppmmin;
end
if(exist('ppmmax','var'))
    ppm_max = ppmmax;
end
if(exist('coilNum','var'))
    coil = coilNum;
end

%center in the dimensions (l-r,a-p,s-i)
center = st.centre;
%dimension size of the MRI
dims = st.vols{1}.dim;
%resolution of MRI (in mm)(not used might be useful)
resolution = (st.bb(2,:) - st.bb(1,:) + 1)/dims(:)';
%3D bounding box of the MRI scan in mm. ie the coordinates where the MRI is
%plotted onto.
bb = st.bb;

range_bool = csi_obj.ppm >= ppm_min & csi_obj.ppm <= ppm_max;

specs = permute(csi_obj.specs, nonzeros([csi_obj.dims.t, csi_obj.dims.x, csi_obj.dims.y, csi_obj.dims.coils, ...
                           csi_obj.dims.averages]));
specs = specs(range_bool, :,:);
ppm = csi_obj.ppm(range_bool);

%get tne yrange for scaling of the y axis to plot the specs 
if(strcmp(type, 'real'))
    yrange = max(max(real(specs(:,:,:,coil)),[],1) - min(real(specs(:,:,:,coil)),[],1), [], 'all');
elseif (strcmp(type, 'imag'))
    yrange = max(max(imag(specs(:,:,:,coil)),[],1) - min(imag(specs(:,:,:,coil)),[],1), [], 'all');
else
    yrange = max(max(abs(specs(:,:,:,coil)),[],1) - min(abs(specs(:,:,:,coil)),[],1), [], 'all');
end
    

%scale factors to fit the spectral dimension at each (x,y) coordinates
scalefactorX=(0.8*csi_obj.deltaX)/(ppm(end) - ppm(1));
scalefactorY=(0.8*csi_obj.deltaY)/yrange;

%set the x dimension and y dimension
xdim = csi_obj.dims.x;
ydim = csi_obj.dims.y;

%reduce the spec scales by scalefactorY
specs = specs(:,:,:,coil).*scalefactorY;

% %reverse the y dimension so plotting anterior = +y and posterior = -y
% tempSpec.specs = flip(tempSpec.specs, ydim);

%set current axis object to be of the axial image. (st.vols{1}.ax{2} and
%st.vols{1}.ax{3} are objects for the sagital and coronal)
axes(st.vols{1}.ax{1}.ax)
%delete previous MRSI plot on MRI

%don't plot unless in the correct z position
%(needs to be modified if 3D MRSI is to be done) 
h = findobj(gca,'Tag','csi_plot');
delete(h);
specs = flip(specs,1);
if(center(3) - csi_obj.deltaZ/2 < csi_obj.imageOrigin(3) && csi_obj.imageOrigin(3) < center(3) + csi_obj.deltaZ/2)
    
    
    %hold axis
        hold on;
        for x = 1:size(csi_obj.specs,xdim)
            for y = 1:size(csi_obj.specs,ydim)
                %first scale the ppm scale so that the range is correct;
                ppm_plot = (ppm-ppm(1))*scalefactorX;
                %now shift it to the correct x-position;
                ppm_plot = ppm_plot + (x-1)*csi_obj.deltaX - (0.8*csi_obj.deltaX)/2 + csi_obj.xCoordinates(1) - bb(1,1);
                %check plottype
                
                %plot whichever plot type is chosen
                switch(type)
                    case 'real'
                        specs = real(specs);
                    case 'imaginary'
                        specs = imag(specs);
                    case 'magnitude'
                        specs = abs(specs);
                    otherwise
                        error("please enter a valid plot_type");
                end
                %Now start plotting
                p = plot(gca, ppm_plot, specs(:,x,y) + (y-1)*csi_obj.deltaY + csi_obj.yCoordinates(1)+ -bb(1,2));
                %set the label of the line object
                set(p, 'Tag', 'csi_plot', 'HitTes', 'off');
            end
        end
    
end
%stop holding
hold off;


end
