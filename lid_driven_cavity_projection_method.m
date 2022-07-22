function [x,y,uf,vf,Ut,pf] = lid_driven_cavity_projection_method()
    clc
    n=256;                  %no of grid points
    dt=0.001;               %time step
    Re=1000;                 %Reynolds no
    tf=60;                  %Final Time
    nx=n; ny=n;
    dx=1/(nx-1); dy=1/(ny-1);
    x=0:dx:1; y=0:dy:1;
    p=zeros(nx+1,ny+1); pn=zeros(nx+1,ny+1); 
    f=zeros(nx+1,ny+1); a=zeros(n+1,n+1);
    u=zeros(nx+1,ny); us=zeros(nx+1,ny);
    v=zeros(nx,ny+1); vs=zeros(nx,ny+1);
    uf=zeros(nx,ny); vf=zeros(nx,ny); pf=zeros(nx,ny); Ut=zeros(nx,ny);
    %specifying initial conditions for u,v, and p
    for i=1:1:n+1
        for j=1:1:n
            u(i,j)=0;
            if i==n+1
                u(i,j)=2;
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
    %constructing 'a' matrix
    for i=1:1:n+1
        for j=1:1:n+1
            if j==i-1||j==i+1   
                a(i,j)=1;
            elseif j==i
                if i==1||i==n+1
                    a(i,j)=-1;
                else
                    a(i,j)=-2;
                end
            end
        end
    end
    [pa,da]=eig(a);
    inpa=inv(pa);
    for k=1:1:tf/dt
        if rem(k,1000)==0
            disp(k*dt) %displays current time
        end
        %x momentum
        for i=2:1:n
            for j=2:1:n-1
                du_2=(u(i,j+1)^2-u(i,j-1)^2)/(2*dx);
                duv=((u(i,j)+u(i+1,j))*(v(i,j)+v(i,j+1))-(u(i,j)+u(i-1,j))*(v(i-1,j)+v(i-1,j+1)))*0.25/dy;
                d2ux=(u(i,j+1)-2*u(i,j)+u(i,j-1))/(dx*dx);
                d2uy=(u(i+1,j)-2*u(i,j)+u(i-1,j))/(dy*dy);
                rhsu=-du_2-duv+(d2ux+d2uy)/Re;
                us(i,j)=u(i,j)+dt*rhsu;
            end
        end
        %x momentum BC
        for j=2:1:n-1
            us(n+1,j)=2-us(n,j);
            us(1,j)=-us(2,j);
        end
        %y momentum
        for i=2:1:n-1
            for j=2:1:n
                dv_2=(v(i+1,j)^2-v((i-1),j)^2)/(2*dy);
                dvu=((v(i,j)+v(i,j+1))*(u(i,j)+u(i+1,j))-(v(i,j-1)+v(i,j))*(u(i+1,j-1)+u(i,j-1)))*0.25/dx;
                d2vx=(v(i,j+1)-2*v(i,j)+v(i,j-1))/(dx*dx);
                d2vy=(v(i+1,j)-2*v(i,j)+v(i-1,j))/(dy*dy);
                rhsv=-dv_2-dvu+(d2vx+d2vy)/Re;
                vs(i,j)=v(i,j)+dt*rhsv;
            end
        end
        %y momentum BC
        for i=2:1:n
            vs(i,1)=-vs(i,2);
            vs(i,n+1)=-vs(i,n);
        end
        %solving pressure poissions eq
        %calc f
        for i=2:1:n
            for j=2:1:n
                f(i,j)=dx*dy*((u(i,j)-u(i,j-1))/(dx)+(v(i,j)-v(i-1,j))/(dy))/(dt);
            end
        end
        fb=inpa*f*pa;
        %cal pressure
        for i=1:1:n+1
            for j=1:1:n+1
                pn(i,j)=fb(i,j)/(da(i,i)+da(j,j));
            end
        end
        pn=pa*pn*inpa;
        %prssure BC
        for j=2:1:n
            pn(1,j)=pn(2,j);
            pn(n+1,j)=pn(n,j);
            pn(j,1)=pn(j,2);
            pn(j,n+1)=pn(j,n);
        end
        %correcting velocities
        for i=2:1:n
            for j=2:1:n-1
                u(i,j)=us(i,j)-dt*(pn(i,j+1)-pn(i,j))/(dx);
            end
        end
        for j=2:1:n-1
            u(n+1,j)=2-u(n,j);
            u(1,j)=-u(2,j);
        end
        for i=2:1:n-1
            for j=2:1:n
                v(i,j)=vs(i,j)-dt*(pn(i+1,j)-pn(i,j))/(dy);
            end
        end
        for i=2:1:n
            v(i,1)=-v(i,2);
            v(i,n+1)=-v(i,n);
        end
        p=pn;
    end
    %cal values on grid points
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
    Ut=sqrt(uf.^2+vf.^2);
    %output
    figure(1);
    contour(x,y,pf,100);
    figure(2);
    quiver(x,y,uf./Ut,vf./Ut,.75);
    figure(3);
    streamslice(x,y,uf,vf);
end