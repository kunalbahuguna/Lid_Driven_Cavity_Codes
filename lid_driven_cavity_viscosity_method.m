function [x,y,pf,uf,vf,phi] = lid_driven_cavity_viscosity_method()
    clc
    n=256;                  %no of grid points
    dt=0.001;               %time step
    Re=1000;                 %Reynolds no
    err_tol=0.0001;     %Error Tolerance
    nx=n; ny=n;
    dx=1/(nx-1); dy=1/(ny-1);
    x=0:dx:1; y=0:dy:1;
    p=zeros(nx+1,ny+1); pn=zeros(nx+1,ny+1); pf=zeros(nx,ny);
    e=zeros(nx+1,ny+1); tot_e=1; c=0; delta=2.5;
    u=zeros(nx+1,ny); us=zeros(nx+1,ny);
    v=zeros(nx,ny+1); vs=zeros(nx,ny+1);
    uf=zeros(nx,ny); vf=zeros(nx,ny);
    f=zeros(nx-2,ny-2); fb=zeros(nx-2,ny-2); phi=zeros(nx-2,ny-2); phib=zeros(nx-2,ny-2);
    for i=1:1:n+1
        for j=1:1:n
            u(i,j)=0;
            if i==n+1
                u(i,j)=1;
            elseif i==n
                u(i,j)=1;
            end
        end
    end
    for i=1:1:n
        for j=1:1:n+1
            v(i,j)=0;
        end
    end
    for i=1:1:n+1
        for j=1:1:n+1
            p(i,j)=1;
        end
    end
    while tot_e>err_tol
        c=c+1;
        %x momentum
        for i=2:1:n
            for j=2:1:n-1
                du_2=(u(i,j+1)^2-u(i,j-1)^2)/(2*dx);
                duv=((u(i,j)+u(i+1,j))*(v(i,j)+v(i,j+1))-(u(i,j)+u(i-1,j))*(v(i-1,j)+v(i-1,j+1)))*0.25/dy;
                d2ux=(u(i,j+1)-2*u(i,j)+u(i,j-1))/(dx*dx);
                d2uy=(u(i+1,j)-2*u(i,j)+u(i-1,j))/(dy*dy);
                p_x=(p(i,j+1)-p(i,j))/dx;
                rhsu=-du_2-duv-p_x+(d2ux+d2uy)/Re;
                us(i,j)=u(i,j)+dt*rhsu;
            end
        end
        %x momentum BC
        for j=2:1:n-1
            us(n+1,j)=2-us(n,j);
            us(1,j)=-us(2,j);
        end
        us;
        %y momentum
        for i=2:1:n-1
            for j=2:1:n
                dv_2=(v(i+1,j)^2-v((i-1),j)^2)/(2*dy);
                dvu=((v(i,j)+v(i,j+1))*(u(i,j)+u(i+1,j))-(v(i,j-1)+v(i,j))*(u(i+1,j-1)+u(i,j-1)))*0.25/dx;
                d2vx=(v(i,j+1)-2*v(i,j)+v(i,j-1))/(dx*dx);
                d2vy=(v(i+1,j)-2*v(i,j)+v(i-1,j))/(dy*dy);
                p_y=(p(i+1,j)-p(i,j))/dy;
                rhsv=-dv_2-dvu-p_y+(d2vx+d2vy)/Re;
                vs(i,j)=v(i,j)+dt*rhsv;
            end
        end
        %y momentum BC
        for i=2:1:n
            vs(i,1)=-vs(i,2);
            vs(i,n+1)=-vs(i,n);
        end
        vs;
        %cal pressure
        for i=2:1:n
            for j=2:1:n
                pn(i,j)=p(i,j)-delta*dt*(  (us(i,j)-us(i,j-1))/dx+(vs(i,j)-vs(i-1,j))/dy  );
            end
        end
        %prssure BC
        for j=2:1:n
            pn(1,j)=pn(2,j);
            pn(n+1,j)=pn(n,j);
            pn(j,1)=pn(j,2);
            pn(j,n+1)=pn(j,n);
        end
        pn;
        tot_e=0;
        for i=2:1:n
            for j=2:1:n
                e(i,j)=(us(i,j)-us(i,j-1))/dx+(vs(i,j)-vs(i-1,j))/dy;
                tot_e=tot_e+abs(e(i,j));
            end
        end
        u=us; v=vs; p=pn; 
        if rem(c,1000)==0
            c
            tot_e
        end
    end
    %cal on grid
    for i=1:1:n-1
        for j=1:1:n-1
            uf(i,j)=0.5*(u(i,j)+u(i+1,j));
        end
    end
    for i=1:1:n-1
        for j=1:1:n-1
            vf(i,j)=0.5*(v(i,j)+v(i,j+1));
        end
    end
    for i=1:1:n-1
        for j=1:1:n-1
            pf(i,j)=.25*(p(i,j)+p(i+1,j)+p(i,j+1)+p(i+1,j+1));
        end
    end
    %stream func
    for i=1:1:n-2
        for j=1:1:n-2
            if j==i-1||j==i+1
                a(i,j)=1;
            elseif j==i
                a(i,j)=-2;
            else
                a(i,j)=0;
            end
        end
    end
    [pa,da]=eig(a);
    inpa=inv(pa);
    for i=1:1:n-2
        for j=1:1:n-2
            f(i,j)=dx*dy*((uf(i+2,j+1)-uf(i,j+1))/(2*dy)-(vf(i+1,j)+vf(i+1,j+2))/(2*dx));
        end
    end
    fb=inpa*f*pa;
    for i=1:1:n-2
        for j=1:1:n-2
            phib(i,j)=fb(i,j)/(da(i,i)+da(j,j));
        end
    end
    phi=pa*phib*inpa;
    Ut=sqrt(uf.^2+vf.^2);
    %output
    figure(1);
    contour(x,y,pf,100);
    figure(2);
    quiver(x,y,uf./Ut,vf./Ut,.75);
    figure(3);
    streamslice(x,y,uf,vf);
end

