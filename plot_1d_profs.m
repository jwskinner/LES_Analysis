function plot_1d_profs(params, xlabels, legendLabels) 
    %% plots the vertical 1d profiles in format (Z, params, data...)

    scrsz = get(0,'ScreenSize');
    figure('Position', [0 0 scrsz(3)*0.9 scrsz(3)*0.9]);

    Z = params{size(params,2)};

    xlimits = zeros(8, 2);
    xlimits(1,:) = [-2, 2];
    xlimits(2,:) = [280, 420];
    xlimits(3,:) = [0, 0.2];
    xlimits(4,:) = [0, 10^-4];
    xlimits(5,:) = [0, 0.5];
    xlimits(6,:) = [0, 1];
    xlimits(7,:) = [0, 10^-15];
    xlimits(8,:) = [-10^-14, 10^-14];

    for i = 1:8
        subplot(2,4,i);
        P1 = plot(params{i}(1,:), Z/10^3, 'LineWidth', 2); hold on; grid on;
        P2 = plot(params{i}(2,:), Z/10^3, 'LineWidth', 2);
        legend([P1, P2], legendLabels, 'Location', 'best')
        xlabel(xlabels{i})
        ylabel('Height [km]')
        xlim([xlimits(i,1), xlimits(i,2)])
    end

    exportgraphics(gcf,strcat('./plots/profiles/', 'profiles', '.gif'),'Resolution',150, 'Append',true)

end
