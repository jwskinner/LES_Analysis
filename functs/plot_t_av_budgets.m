function plot_t_av_budgets(budget_array, scaling, Z, i_start, i_end, n_budgets, budgets, ptitle)
    
    % This function plots the budget terms as computed by the time avering
    % budget plots (plot_tavg_XX_budgets.m)
    
    % Example usage
    % uv_t = ...; % Provide the uv_t data (replace ... with your data)
    % Z = ...; % Provide the Z data (replace ... with your data)
    % i_start = ...; % Provide the starting index for the loop (replace ... with your value)
    % i_end = ...; % Provide the ending index for the loop (replace ... with your value)
    % n_budgets = ...; % Provide the number of budgets (replace ... with your value)
    % budgets = {...}; % Provide the budgets cell array (replace ... with your budgets)

    % Some parameters to vary 
    xlims = [-5e-4 5e-4]; 
    ylims = [0 18]; 

    
    % PLOT ALL THE BUDGET TERMS
    figure('Renderer', 'painters', 'Position', [10 10 900 600])            % makes paper format figure
    lw = 1.5;
    colors = {'#0072BD', '#A2142F', '#EDB120', '#77AC30', '#4DBEEE', '#800080', '#000000'}; 
    style = repmat({'-'}, 1, n_budgets-1);                                 % Generate n-1 occurrences of '-' for the budget line styles 
    style = [style, '--'];                                                 % The residuals are dashed '--'                  
    alpha = 1.0; 

    % DIAGNOSTIC PLOT OF TIME-VARYING BUDGETS
    rgbC = [];
    for i = 1:n_budgets 
        rgb = sscanf(colors{i}(2:end), '%2x')/255;
        rgbC = [rgbC, rgb]; 
    end

    subplot(1, 1, 1)
    for i = 1:(i_end-i_start)
        for j = 1:n_budgets
            plot(squeeze(budget_array(i, j, :)), Z/10^3, style{j}, ...
                'LineWidth', lw, 'Color', [rgbC(:,j)', alpha]); hold on;
        end
    end
    ylabel('$z$ [km]', 'LineWidth', 1.5, 'FontSize', 15, 'Interpreter', ...
        'Latex')
    legend(budgets);
    title([ptitle,' variance budget'], 'FontSize', 15, ...
        'Interpreter', 'Latex'); hold off;
    xlim(xlims)
    ylim(ylims)

    % MAIN PLOT OF TIME AVERAGED BUDGETS
    figure('Renderer', 'painters', 'Position', [10 10 900 600])            % makes paper format figure
    subplot(1, 1, 1)
    budget_array_mean = squeeze(mean(budget_array, 1));                    % Calculate mean budgets over the time dimension
    for i = 1:n_budgets
        plot(budget_array_mean(i, :), Z/10^3, style{i}, 'LineWidth', lw,...
            'Color', colors{i}); hold on;
    end 
    ylabel('$z$ [km]', 'LineWidth', 1.5, 'FontSize', 15, 'Interpreter', ...
        'Latex')
    legend(budgets);
    title(['$t$ average ', ptitle, ' variance budget'], 'FontSize', 15, ...
        'Interpreter', 'Latex'); hold off;
    xlim(xlims)
    ylim(ylims)

    % SUBPLOT VERSION OF TIME AVERAGED BUDGETS
    figure('Renderer', 'painters', 'Position', [10 10 900 600])            % makes paper format figure
    for i = 1:n_budgets
        subplot(1, n_budgets, i);                                          % Create a new subplot for each line
        plot(budget_array_mean(i, :), Z/10^3, 'k', 'LineWidth', lw); 
        ylabel('$z$ [km]', 'LineWidth', 1.5, 'FontSize', 15, ...
            'Interpreter', 'Latex')
        title(['Line ' num2str(i)], 'FontSize', 12);                       % Add title for each subplot (line)
        xlim(xlims)
        ylim(ylims)
        title(budgets{i}, 'FontSize', 12);                                 % Add title for the first subplot
        sgtitle([ptitle], 'FontSize', 15, 'Interpreter', ...
            'Latex');
    end

     % SUBPLOT VERSION OF TIME AVERAGED BUDGETS WITH SCALING
    figure('Renderer', 'painters', 'Position', [10 10 900 600])            % makes paper format figure
    for i = 1:n_budgets
        subplot(1, n_budgets, i);                                          % Create a new subplot for each line
        plot(budget_array_mean(i, :)/scaling.x, scaling.y, 'k', 'LineWidth', lw); 
        ylabel('$z/z^*$', 'LineWidth', 1.5, 'FontSize', 15, ...
            'Interpreter', 'Latex')
        xlabel('$\rho_* w_*^3 / z^*$', 'LineWidth', 1.5, 'FontSize', 15, ...
            'Interpreter', 'Latex')
        title(['Line ' num2str(i)], 'FontSize', 12);                       % Add title for each subplot (line)
        xlim([-1, 1])
        ylim([0, 1.5])
        title(budgets{i}, 'FontSize', 12);                                 % Add title for the first subplot
        sgtitle([ptitle], 'FontSize', 15, 'Interpreter', ...
            'Latex');
    end
end