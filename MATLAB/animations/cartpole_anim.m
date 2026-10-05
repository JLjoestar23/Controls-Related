function cartpole_anim(X, t, tspan, system_params)

    l  = system_params.l;
    mc = system_params.mc;
    mp = system_params.mp;

    cart_w = 0.3;   % cart width  (m)
    cart_h = 0.15;  % cart height (m)

    % cartesian position of pendulum bob
    x_cart =  X(:,1);
    x_bob  =  X(:,1) + l*sin(X(:,2));
    y_bob  = -l*cos(X(:,2));       % 0 = pivot height

    % axis limits
    x_max = max(x_cart)*1.1 + l;
    x_min = min(x_cart)*1.1 - l;

    % ── figure ────────────────────────────────────────────────────────────
    fig = figure('Name','Cart-Pole');

    % left panel: animation
    ax1 = subplot(1,2,1,'Parent',fig);
    set(ax1,'FontSize',10);
    axis(ax1,'equal'); grid(ax1,'on');
    xlim(ax1,[x_min  x_max]);
    ylim(ax1,[-(l*1.3)   l*1.3]);
    xlabel(ax1,'x (m)'); ylabel(ax1,'y (m)');
    title(ax1,'Cart-Pole','FontSize',12,'FontWeight','bold');
    hold(ax1,'on');

    % ground rail
    plot(ax1,[x_min x_max],[0 0],'k-','LineWidth',1.5);

    % trail of pendulum bob
    trail_len = 200;
    h_trail = plot(ax1, NaN, NaN, '-', 'Color',[0.2 0.6 1 0.6], 'LineWidth',1.2);

    % cart (drawn as a filled rectangle via patch)
    cx = x_cart(1);
    cart_vx = cx + [-1 1 1 -1]*cart_w/2;
    cart_vy = [0 0 1 1]*cart_h - cart_h/2;
    h_cart  = patch(ax1, cart_vx, cart_vy, [0.4 0.4 0.4], 'EdgeColor','k','LineWidth',1.2);

    % rod
    h_rod = plot(ax1,[cx x_bob(1)],[0 y_bob(1)],'-','Color',[0.75 0.75 0.75],'LineWidth',2.5);

    % pivot dot (moves with cart)
    h_pivot = plot(ax1, cx, 0, 'o','MarkerSize',7, ...
                   'MarkerFaceColor','k','MarkerEdgeColor','w');
    % pendulum bob
    % h_bob = plot(ax1, 0, 0, 'o','MarkerSize',18, ...
    %              'MarkerFaceColor',[0.2 0.6 1],'MarkerEdgeColor','k','LineWidth',1);
    radius = 0.05;
    h_bob = rectangle(ax1, 'Position', [x_bob(1)-radius, y_bob(1)-radius, ...
                                     2*radius, 2*radius], ...
                    'Curvature', [1 1], ...
                    'FaceColor', [0.2 0.6 1], ...
                    'EdgeColor', 'k', ...
                    'LineWidth', 1);

    % time label
    h_time = text(ax1, x_min*0.93, l*1.15, 't = 0.00 s','Color','k','FontSize',13);

    % right panel: phase portrait (theta vs dtheta)
    ax2 = subplot(1,2,2,'Parent',fig);
    set(ax2,'FontSize',10);
    axis(ax2,'square'); grid(ax2,'on');
    th_all = wrapToPi(X(:,2));
    xlim(ax2,[min(th_all)*1.2, max(th_all)*1.2]);
    ylim(ax2,[min(X(:,4))*1.2, max(X(:,4))*1.2]);
    xlabel(ax2,'$\theta$ (rad)','Interpreter','latex');
    ylabel(ax2,'$\dot{\theta}$ (rad/s)','Interpreter','latex');
    title(ax2,'Pendulum Phase Portrait','FontSize',12,'FontWeight','bold');
    hold(ax2,'on');
    pp = animatedline(ax2,'Color',[0.2 0.6 1],'LineWidth',1.5);

    % ── interpolate to fixed frame rate ───────────────────────────────────
    frame_rate = 48;
    t_frames   = tspan(1):1/frame_rate:tspan(2);

    xf      = interp1(t, X(:,1), t_frames);
    thetaf  = interp1(t, X(:,2), t_frames);
    dthetaf = interp1(t, X(:,4), t_frames);
    xbf     = interp1(t, x_bob,  t_frames);
    ybf     = interp1(t, y_bob,  t_frames);

    trail_x = nan(1, trail_len);
    trail_y = nan(1, trail_len);

    pause(1);

    % ── animation loop ─────────────────────────────────────────────────────
    tic;
    for k = 1:length(t_frames)

        % real-time pacing
        elapsed = toc;
        if elapsed < t_frames(k)
            pause(t_frames(k) - elapsed);
        end

        % trail (circular buffer)
        idx = mod(k-1, trail_len) + 1;
        trail_x(idx) = xbf(k);
        trail_y(idx) = ybf(k);
        if k <= trail_len
            tx = trail_x(1:k);
            ty = trail_y(1:k);
        else
            order = [idx+1:trail_len, 1:idx];
            tx = trail_x(order);
            ty = trail_y(order);
        end

        % cart vertices
        cvx = xf(k) + [-1 1 1 -1]*cart_w/2;
        cvy = [0 0 1 1]*cart_h - cart_h/2;

        set(h_trail,  'XData', tx,               'YData', ty);
        set(h_rod,    'XData', [xf(k), xbf(k)],  'YData', [0,    ybf(k)]);
        % set(h_bob,    'XData', xbf(k),            'YData', ybf(k));
        
        h_bob.Position(1) = xbf(k) - radius;
        h_bob.Position(2) = ybf(k) - radius;

        set(h_cart,   'XData', cvx,               'YData', cvy);
        set(h_pivot,  'XData', xf(k),             'YData', 0);
        set(h_time,   'String', sprintf('t = %.2f s', t_frames(k)));

        addpoints(pp, thetaf(k), dthetaf(k));

        drawnow limitrate;
        if ~ishandle(fig), break; end

    end

end