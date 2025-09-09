%%%
   %计算二维矩形平面NS方程 
   %有限差分法
%%%
    % 参数设置
    nx = 101;    
    ny = 101;
    nt = 1000;    %迭代次数
    nit = 1000;   %求解压强迭代次数
    rho = 1;
    nu = 0.001;
    dt = 0.001;  %时间步
    
    dx = 2/(nx-1);
    dy = 2/(ny-1);
    
    % 初始化物理场
    u = zeros(ny,nx);
    v = zeros(ny,nx);
    p = zeros(ny,nx);
    b = zeros(ny,nx);
    
    % 网格坐标划分
    x = linspace(0,2,nx);
    y = linspace(0,2,ny);
    [X,Y] = meshgrid(x,y);
    
    % 主时间循环  1-nt
    for n = 1:nt
        un = u;  %存储上一步的X方向速度场
        vn = v;  %存储上一步的Y方向速度场
        
        % 计算泊松方程源项 b
        b(2:end-1,2:end-1) = rho*(1/dt*...
            ((u(2:end-1,3:end)-u(2:end-1,1:end-2))/(2*dx) + ...
             (v(3:end,2:end-1)-v(1:end-2,2:end-1))/(2*dy)) - ...
            ((u(2:end-1,3:end)-u(2:end-1,1:end-2))/(2*dx)).^2 - ...
            2*((u(3:end,2:end-1)-u(1:end-2,2:end-1))/(2*dy).*...
               (v(2:end-1,3:end)-v(2:end-1,1:end-2))/(2*dx)) - ...
            ((v(3:end,2:end-1)-v(1:end-2,2:end-1))/(2*dy)).^2);
        
        % 压强泊松方程
        pn = p;
        for q = 1:nit
            pn = p;   %存储上一步的压强场
            p(2:end-1,2:end-1) = ((pn(2:end-1,3:end)+pn(2:end-1,1:end-2))*dy^2 + ...
                                  (pn(3:end,2:end-1)+pn(1:end-2,2:end-1))*dx^2) / ...
                                 (2*(dx^2+dy^2)) - ...
                                 dx^2*dy^2/(2*(dx^2+dy^2))*b(2:end-1,2:end-1);
            
            % 边界条件
            p(:,end) = p(:,end-1);   % dp/dx = 0 at x = 2
            p(1,:) = p(2,:);         % dp/dy = 0 at y = 0
            p(:,1) = p(:,2);         % dp/dx = 0 at x = 0
            p(end,:) = 0;            % p = 0 at y = 2
        end
        
        % 更新速度场  用当前的压强场 求解下一步的速度场
        u(2:end-1,2:end-1) = un(2:end-1,2:end-1) - ...
            un(2:end-1,2:end-1).*dt/dx.*(un(2:end-1,2:end-1)-un(2:end-1,1:end-2)) - ...
            vn(2:end-1,2:end-1).*dt/dy.*(un(2:end-1,2:end-1)-un(1:end-2,2:end-1)) - ...
            dt/(2*rho*dx)*(p(2:end-1,3:end)-p(2:end-1,1:end-2)) + ...
            nu*(dt/dx^2*(un(2:end-1,3:end)-2*un(2:end-1,2:end-1)+un(2:end-1,1:end-2)) + ...
                dt/dy^2*(un(3:end,2:end-1)-2*un(2:end-1,2:end-1)+un(1:end-2,2:end-1)));
            
        v(2:end-1,2:end-1) = vn(2:end-1,2:end-1) - ...
            un(2:end-1,2:end-1).*dt/dx.*(vn(2:end-1,2:end-1)-vn(2:end-1,1:end-2)) - ...
            vn(2:end-1,2:end-1).*dt/dy.*(vn(2:end-1,2:end-1)-vn(1:end-2,2:end-1)) - ...
            dt/(2*rho*dy)*(p(3:end,2:end-1)-p(1:end-2,2:end-1)) + ...
            nu*(dt/dx^2*(vn(2:end-1,3:end)-2*vn(2:end-1,2:end-1)+vn(2:end-1,1:end-2)) + ...
                dt/dy^2*(vn(3:end,2:end-1)-2*vn(2:end-1,2:end-1)+vn(1:end-2,2:end-1)));
        
        % 速度边界条件
        u(1,:) = 0;      % u = 0 at y = 0
        u(:,1) = 0;      % u = 0 at x = 0
        u(:,end) = 0;    % u = 0 at x = 2
        u(end,:) = 0.5;    % u = 1 at y = 2 
        v(1,:) = 0;      % v = 0 at y = 0
        v(end,:) = 0;    % v = 0 at y = 2
        v(:,1) = 0;      % v = 0 at x = 0
        v(:,end) =0 ;    % v = 0 at x = 2
    end
    
    % 可视化处理
    figure;
    contourf(X,Y,p,20,'LineColor','none');
    colorbar;
    hold on;
    quiver(X(1:2:end,1:2:end),Y(1:2:end,1:2:end),u(1:2:end,1:2:end),v(1:2:end,1:2:end));
    xlabel('X');
    ylabel('Y');
    title('Cavity Flow');
    hold off;

