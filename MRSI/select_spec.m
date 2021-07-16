%select_spec.m 
%
% callback function used for uicontrol popup menu for choosing to select
% MRSI specs to plot. Plots the closest clicked MRSI spec in spm_image GUI.
% must be used with sim_CSIoverlayMRI and spm_image. See
% https://www.mathworks.com/help/matlab/creating_plots/callback-definition.html
% for more info on callback funcitons.
%
% USAGE: can't be used by itself, must be called back
% 
% INPUT: 
% src       = modified uicontrol object
%

function select_spec(src,~) 
    %get global MRSI object
    global csi_obj;
    global ppm_min;
    global ppm_max;
    global type;
    
    %only start if the users clicks with the right mouse button
    if ~strcmpi(get(gcbf,'SelectionType'),'alt')
        %get where the user clicked on the plot
        mousePos=get(src,'CurrentPoint');
        %get x and y coordinates
        mousePos=mousePos(1:2:3);
        %get all the MRSI line plots
        children = findobj(src,'Tag','csi_plot');
        %flip the array so the spec plotted first in overlay is first in
        %the array
        children = flipud(children);
        %temp values
        min = realmax;
        closest = -1;
        for i = 1: numel(children)
            %get the distance from the specs plot to the clicked position
            cur = children(i);
            position = [mean(cur.XData) mean(cur.YData)];
            distance = sqrt(sum((mousePos-position).^2));
            
            %change the temp values if the distance is smaller than min
            if(distance < min)
                min = distance;
                closest = i;
            end
        end
        %display the linear index of the closest MRSI spec
        disp(closest)
        %get the corresponding x,y indecies from linear indexing
        [row, col] = ind2sub([csi_obj.sz(csi_obj.dims.y), csi_obj.sz(csi_obj.dims.x)], closest);
        %invert the y coordinates
        row = csi_obj.sz(csi_obj.dims.y) - row + 1;
        
        %delete the current MRSI shown plot if there is one
        ax = findobj(gcf, 'Tag', 'csi_axis');
        if(numel(ax.Children) ~= 0)
            delete(ax.Children)
        end
        %hold the axis
        hold(ax, 'on')
        %plot the selected MRSI spec.4
        range_bool = csi_obj.ppm > ppm_min & csi_obj.ppm < ppm_max;
        ppm = csi_obj.ppm(range_bool);
        set(ax,'XDir','reverse')
        xlim(ax, [ppm_min ppm_max]);
        
        switch(type)
            case 'real'
                specs = real(csi_obj.specs(range_bool, col, row));
            case 'imaginary'
                specs = imag(csi_obj.specs(range_bool, col, row));
            case 'magnitude'
                specs = abs(csi_obj.specs(range_bool, col, row));
            otherwise
                error("please enter a valid plot_type");
        end
        
        plot(ax, ppm, specs)
        hold(ax, 'off')
        
    end 

end