function plot_1d_profs(params, xlabels, legendLabels, time) 
    %% plots the vertical 1d profiles in format (Z, params, data...)

    scrsz = get(0,'ScreenSize');
    figure('Position', [0 0 scrsz(3)*0.9 scrsz(3)*0.9]);

    Z = params{size(params,2)};

    xlimits = zeros(8, 2);
    xlimits(1,:) = [-2, 2];
    xlimits(2,:) = [290, 320];
    xlimits(3,:) = [0.005, 0.02];
    xlimits(4,:) = [0, 0.2*10^-4];
    xlimits(5,:) = [0, 0.2];
    xlimits(6,:) = [0, 1];
    xlimits(7,:) = [0, 0.2*10^-15];
    xlimits(8,:) = [-10^-14, 10^-14];


    ylimits = [0,4];

    for i = 1:8
        subplot(2,4,i);
        for j = 1:size(params{i}, 1) % Loops over lines to plot
            P(i) = plot(params{i}(j,:), Z/10^3, 'LineWidth', 2); ...
                hold on; grid on;
        end
        legend(legendLabels, 'Location', 'northeast')
        xlabel(xlabels{i})
        ylabel('Height [km]')
        xlim([xlimits(i,1), xlimits(i,2)])
        ylim(ylimits)
    end
    sgtitle(strcat(num2str(time), ' hours'))
    exportgraphics(gcf,strcat('./plots/profiles/', 'profiles', '.gif'),...
        'Resolution',150, 'Append',true)

end
