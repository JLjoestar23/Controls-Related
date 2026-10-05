function double_pendulum_anim(X, t, tspan, system_params)

    l1 = system_params.l1;
    l2 = system_params.l2;

    % convert angles to cartesian coordinates
    x1 =  l1*sin(X(:,1));
    y1 = -l1*cos(X(:,1));
    x2 =  x1 + l2*sin(X(:,2));
    y2 =  y1 - l2*cos(X(:,2));
    
    % figure setup
    fig = figure('Name','Double Pendulum');
    
    % left panel: animation
    ax1 = subplot(1,2,1,'Parent',fig);
    set(ax1,'FontSize',10);
    axis(ax1,'equal'); grid(ax1, 'on');
    lim = (l1+l2)*1.25;
    xlim(ax1,[-lim lim]); ylim(ax1,[-lim lim]);
    xlabel(ax1,'x (m)'); ylabel(ax1,'y (m)');
    title(ax1,'Double Pendulum','FontSize',12,'FontWeight','bold');
    hold(ax1,'on');
    
    % trail line
    trail_len = 200;
    h_trail = plot(ax1, NaN, NaN, '-', 'Color',[0.2 0.6 1 0.6], 'LineWidth',1.2);
    
    % rods
    h_rod1  = plot(ax1,[0 0],[0 0],'-','Color',[0.75 0.75 0.75],'LineWidth',2.5);
    h_rod2  = plot(ax1,[0 0],[0 0],'-','Color',[0.75 0.75 0.75],'LineWidth',2.5);
    
    % pivot
    plot(ax1, 0, 0, 'o','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','w');
    
    % masses
    h_m1 = plot(ax1,0,0,'o','MarkerSize',18,'MarkerFaceColor',[0.2 0.6 1], ...
                'MarkerEdgeColor','k','LineWidth',1);
    h_m2 = plot(ax1,0,0,'o','MarkerSize',18,'MarkerFaceColor',[1 0.35 0.18], ...
                'MarkerEdgeColor','k','LineWidth',1);
    
    % time label
    h_time = text(ax1,-lim*0.95, lim*0.88,'t = 0.00 s','Color','k','FontSize',15);
    
    % right panel: angle time series (live)
    % ax2 = subplot(1,2,2,'Parent',fig);
    % set(ax2,'FontSize',10); 
    % axis(ax2,'square'); grid(ax2,'on');
    % xlim(ax2,[0 tspan(2)]); ylim(ax2,[min([X(:,1); X(:,2)])*1.2 max([X(:,1); X(:,2)])*1.2]);
    % xlabel(ax2,'Time (s)'); ylabel(ax2,'Angle (rad)');
    % title(ax2,'Joint Angles','FontSize',12,'FontWeight','bold');
    % hold(ax2,'on');
    % h_p1plot = animatedline(ax2,'Color',[0.2 0.6 1],'Marker', '.','MarkerSize',10);
    % h_p2plot = animatedline(ax2,'Color',[1 0.35 0.18],'Marker', '.','MarkerSize',10);
    % legend(ax2,{'\phi_1','\phi_2'},'Location','northeast','FontSize',10);
    
    % with this:
    ax3 = subplot(1,2,2,'Parent',fig);
    set(ax3,'FontSize',10);
    axis(ax3,'square'); grid(ax3,'on');
    xlim(ax3,[min(X(:,1))*1.2, max(X(:,1))*1.2]);
    ylim(ax3,[min(X(:,2))*1.2, max(X(:,2))*1.2]);
    xlabel(ax3,'\phi_1 (rad)'); ylabel(ax3,'\phi_2 (rad)');
    title(ax3,'Phase Portrait','FontSize',12,'FontWeight','bold');
    hold(ax3,'on');
    pp = animatedline(ax3,'Color',[0.2 0.6 1], 'LineWidth', 2);
    
    % animation loop
    frame_rate = 48;
    dt_frame   = 1/frame_rate;
    t_frames   = tspan(1):dt_frame:tspan(2);
    
    % interpolate to fixed frame times
    phi1f = interp1(t, X(:,1), t_frames);
    phi2f = interp1(t, X(:,2), t_frames);
    x1f   = interp1(t, x1,     t_frames);
    y1f   = interp1(t, y1,     t_frames);
    x2f   = interp1(t, x2,     t_frames);
    y2f   = interp1(t, y2,     t_frames);
    
    trail_x = nan(1, trail_len);
    trail_y = nan(1, trail_len);
    
    tic;
    for k = 1:length(t_frames)
        % keep real-time pacing
        elapsed = toc;
        target  = t_frames(k);
        if elapsed < target
            pause(target - elapsed);
        end
    
        % update trail
        idx = mod(k-1, trail_len) + 1;
        trail_x(idx) = x2f(k);
        trail_y(idx) = y2f(k);
    
        if k <= trail_len
            tx = trail_x(1:k);
            ty = trail_y(1:k);
        else
            order = [idx+1:trail_len, 1:idx];
            tx = trail_x(order);
            ty = trail_y(order);
        end
    
        set(h_trail, 'XData', tx, 'YData', ty);
        set(h_rod1,  'XData', [0, x1f(k)],       'YData', [0, y1f(k)]);
        set(h_rod2,  'XData', [x1f(k), x2f(k)],  'YData', [y1f(k), y2f(k)]);
        set(h_m1,    'XData', x1f(k), 'YData', y1f(k));
        set(h_m2,    'XData', x2f(k), 'YData', y2f(k));
        set(h_time,  'String', sprintf('t = %.2f s', t_frames(k)));
    
        % addpoints(h_p1plot, t_frames(k), phi1f(k));
        % addpoints(h_p2plot, t_frames(k), phi2f(k));
        addpoints(pp, phi1f(k), phi2f(k));
    
        drawnow limitrate;
    
        if ~ishandle(fig), break; end
    end

end