% Jack W. Skinner generic funciton for plotting 2D fields from LES Data.
% Last touched 6/2/2023

function plot_fields(varargin)

    nar = nargin-2;
    nam = varargin{nar + 1};
    params = varargin{nar + 2};

    subplot_row = ceil(nar/5); 
    subplot_col = mod(nar-1,5)+1; 

    % Get the screen size
    scrsz = get(0,'ScreenSize');

    % Set the figure width to be equal to the screen width
    clear figure
    figure('Position', [0 0 scrsz(3)*0.5 scrsz(3)*0.5]);
    
    for i = 1:nar
        data = varargin{i}.Data;

        s=size(data);n=s(1); m=s(2); nm=n*m;

        if isfield(params, 'mean') 
            data = mean(data, 3); % Mean over the vertical 
        end 

        if isfield(params, 'sum') 
            data = sum(data, 3); % Mean over the vertical 
        end 

        if isfield(params, 'absum') 
            data = sum(abs(data), 3); % Mean over the vertical 
        end 

        if isfield(params, 'autoscale')
            params.cmax = max(data(:));
            params.cmin = min(data(:));
        end 
        
        
        size(data)

        x_km = linspace(0, (n*nam.dx)/1000, n);
        y_km = linspace(0, (n*nam.dy)/1000, n);

        subplot_i= subplot(subplot_row, subplot_col, i);

        imagesc(x_km, y_km, data);
        set(gca,'YDir','normal'); % Flip the y-axis
        ylabel('y [km]','LineWidth',1,'FontSize',13);
        xlabel('x [km]','LineWidth',1,'FontSize',13);
        axis image;
        
        % Format titles from NETCDF Descriptions 
        title_str=varargin{i}.Desc;
        title_str=lower(title_str);
        idx=regexp([' ' title_str],'(?<=\s+)\S','start')-1;
        title_str(idx)=upper(title_str(idx));
        title(append(title_str, ', ', nam.txt),'FontWeight','Normal');

        box(subplot_i,'on');
        axis(subplot_i,'tight');
        set(subplot_i,'DataAspectRatio',[1 1 1],'FontSize',13,'Layer','top',...
        'LineWidth',1,'XLimitMethod','tight','YLimitMethod','tight','ZLimitMethod',...
        'tight');

        c = colorbar;
        c.Label.String = append('[', varargin{i}.Units, ']');
        
        if isfield(params, 'cmin') && isfield(params, 'cmax')
            caxis([params.cmin, params.cmax]);
        end
    end

    sgtitle(append(string(params.time), ' hours')); 

    % Save figure to a folder if specified
    if isfield(params, 'save_folder')
        if ~exist(params.save_folder, 'dir')
            mkdir(params.save_folder);
        end
        exportgraphics(gcf, strcat(params.save_folder, nam.txt, '_', params.name, ".png"),'Resolution',150);
    end

    if isfield(params, 'save_movie') 
     exportgraphics(gcf,strcat(params.save_folder, params.save_movie, '_', varargin{i}.Name, '.gif'),'Resolution',150, 'Append',true)
    end 
end
