	program fdtd

	implicit none

    integer nx, nz, center_x, center_z
    parameter (nx = 500, nz = 500)
    integer xind, zind, i, k, n
    integer max_iter
    parameter (max_iter = 10000)
    
    real*8 l(nx,nz), m(nx,nz)
    real*8 rho(nx,nz), b(nx,nz)
    real*8 c, cmax
    real*8 t0, width, freq, omega, pi, lambda
    real*8 dx, dz, dt, tw
    real*8 tau_xx_det(max_iter), tau_zz_det(max_iter), u_det(max_iter)
    
    real*8 tau_xx(nx,nz)       
    real*8 tau_zz(nx,nz)
    real*8 tau_xz(nx-1,nz-1)
    real*8 ux(nx-1,nz)    
    real*8 uz(nx,nz-1)    

    real*8 tau_xx_n(nx,nz)       
    real*8 tau_zz_n(nx,nz)       
    real*8 ux_n(nx-1,nz)    
    real*8 uz_n(nx,nz-1)    

    complex*16 tau_xx_mod(nx,nz)
    complex*16 tau_zz_mod(nx,nz)       
    complex*16 tau_xz_mod(nx-1,nz-1)       
    complex*16 ux_mod(nx-1,nz)    
    complex*16 uz_mod(nx,nz-1)    

    complex*16 COMPLEX_I    
    real*8 tbplt(nx-1,nz-1), div_x, div_z, xdist, zdist, dist
    character*32 filename
   
    tau_xx = 0.0       
    tau_zz = 0.0
    tau_xz = 0.0
    ux = 0.0    
    uz = 0.0    
   
    COMPLEX_I = dcmplx(0.0, 1.0)

    t0 = 500.0       !index time of gaussian pulse peak
    width = 1000.0     !peakyness of gaussian pulse

    pi = 3.1415926
    freq = 1.7e6
    cmax = 1600.0
    omega = 2*pi*freq
    lambda = cmax/freq

    dx = 5e-6
    dz = dx

    !! Initialise the p and u fields and set coefficients + a Courrant compliant dt
        
    dt = dx/(dsqrt(2.0D0)*cmax) !set time step small enough to satisfy courrant condition
    dt = 0.95*dt
    
    write(*,*) 'dt = ', dt

    !definition of the waveguide structure
    center_x = int(nx/2)-1
    center_z = int(nz/2)-1
    
    rho = 1000.0
    l = rho*(1540.0**2)
    m = 0.0       

    do zind = 1, nz
    do xind = 1, nx
     
        zdist = (center_z-zind)
        xdist = (center_x-xind)
        dist = dsqrt(xdist**2+zdist**2)
        
        if(dist<525e-6/dx/2.0) then
        
		l(xind,zind) = 1100.0*(1400.0**2 - 2.0*440.0**2)
		m(xind,zind) = 1100.0*440.0**2       			
		rho(xind,zind) = 1100.0
    
        endif
        
    enddo
    enddo
    
    do zind = 1, nz
    do xind = 1, nx
     
        zdist = (center_z-zind)
        xdist = (center_x-xind)
        dist = dsqrt(xdist**2+zdist**2)
        
        if(dist<250e-6/dx/2.0) then
        
		l(xind,zind) = 1000.0*1540.0**2
		m(xind,zind) = 0.0
		rho(xind,zind) = 1000.0
    
        endif
        
    enddo
    enddo

	b = 1.0/rho 
       
    !! FDTD loop - comment out plot lines for speed
    do n = 0,max_iter
    
    ! Update the ux values:
    do zind=2,nz-1
    do xind=2,nx-1
    
    ux(xind,zind) = ux(xind,zind) &
    + dt*b(xind,zind)*((tau_xx(xind,zind) - tau_xx(xind-1,zind))/dx &
    + (tau_xz(xind,zind) - tau_xz(xind,zind-1))/dz)
    enddo
    enddo

    do zind=2,nz-1
    do xind=2,nx-1

    uz(xind,zind) = uz(xind,zind) &
    + dt*b(xind,zind)*((tau_xz(xind,zind) - tau_xz(xind-1,zind))/dx &
    + (tau_zz(xind,zind) - tau_zz(xind,zind-1))/dz)
 
    enddo
    enddo

    do zind=1,nz-2
    do xind=1,nx-2
    
    tau_xx(xind,zind) = tau_xx(xind,zind) + dt*(l(xind,zind)+2.0*m(xind,zind))*((ux(xind+1,zind) - ux(xind,zind))/dx) + &
    dt*l(xind,zind)*(uz(xind,zind+1) - uz(xind,zind))/dz

    tau_zz(xind,zind) = tau_zz(xind,zind) + dt*(l(xind,zind)+2.0*m(xind,zind))*((uz(xind,zind+1) - uz(xind,zind))/dz) + &
    dt*l(xind,zind)*(ux(xind+1,zind) - ux(xind,zind))/dx

    enddo
    enddo
    
    do zind=1,nz-1
    do xind=1,nx-1
    
    tau_xz(xind,zind) = tau_xz(xind,zind) + dt*m(xind,zind)*((ux(xind,zind+1) - ux(xind,zind))/dz &
    + (uz(xind+1,zind) - uz(xind,zind))/dz)

    enddo
    enddo

    tau_xx(center_x,5) = tau_xx(center_x,5) + dsin(2.0*pi*n*dt*freq) 
    tau_zz(center_x,5) = tau_zz(center_x,5) + dsin(2.0*pi*n*dt*freq) 

        
    !Mur1 ABC 

    c = 1540.0
   

	do zind=1,nz
    do xind=1,nx

	!LEFT WALL:
        tau_xx(1,zind)=tau_xx_n(2,zind)+((c*dt-dx)/(c*dt+dx))*(tau_xx(2,zind)-tau_xx_n(1,zind))
        tau_zz(1,zind)=tau_zz_n(2,zind)+((c*dt-dx)/(c*dt+dx))*(tau_zz(2,zind)-tau_zz_n(1,zind))

        !RIGHT WALL:
        tau_xx(nx-1,zind)=tau_xx_n(nx-1-1,zind)+((c*dt-dx)/(c*dt+dx))*(tau_xx(nx-1-1,zind)-tau_xx_n(nx-1,zind))
        tau_zz(nx-1,zind)=tau_zz_n(nx-1-1,zind)+((c*dt-dx)/(c*dt+dx))*(tau_zz(nx-1-1,zind)-tau_zz_n(nx-1,zind))

        !TOP WALL
        tau_xx(xind,1)=tau_xx_n(xind,2)+((c*dt-dz)/(c*dt+dz))*(tau_xx(xind,2)-tau_xx_n(xind,1))
        tau_zz(xind,1)=tau_zz_n(xind,2)+((c*dt-dz)/(c*dt+dz))*(tau_zz(xind,2)-tau_zz_n(xind,1))

        !BOTTOM WALL
        tau_xx(xind,nz-1)=tau_xx_n(xind,nz-1-1)+((c*dt-dz)/(c*dt+dz))*(tau_xx(xind,nz-1-1)-tau_xx_n(xind,nz-1))
        tau_zz(xind,nz-1)=tau_zz_n(xind,nz-1-1)+((c*dt-dz)/(c*dt+dz))*(tau_zz(xind,nz-1-1)-tau_zz_n(xind,nz-1))       

    enddo
    enddo

    !save 'old' fields for Mur ABC
    ux_n=ux
    uz_n=uz
    tau_xx_n = tau_xx
    tau_zz_n = tau_zz
    
    !symmetry
    ux(center_x:nx-1,1:nz) = 0.0
    uz(center_x:nx,1:nz-1) = 0.0
    tau_xx(center_x:nx,1:nz) = 0.0
    tau_zz(center_x:nx,1:nz) = 0.0  
    tau_xz(center_x:nx-1,1:nz-1) = 0.0      
              
    if ( mod(n,int(0.1e-6/dt)).eq.0) then   
       
        tau_xx_mod = tau_xx_mod + tau_xx*exp(COMPLEX_I*n*dt*2.0*pi*freq)
		tau_zz_mod = tau_zz_mod + tau_zz*exp(COMPLEX_I*n*dt*2.0*pi*freq)
		ux_mod = ux_mod + ux*exp(COMPLEX_I*n*dt*2.0*pi*freq)
		uz_mod = uz_mod + uz*exp(COMPLEX_I*n*dt*2.0*pi*freq)
        
         
    endif

	if ( mod(n,100).eq.0) then   

        write(*,*) 'step #', n, 'time ', dble(n)*dt/1e-6, 'micro sec'
                
        tbplt = dsqrt(abs(ux_mod(1:nx-1,1:nz-1))**2+abs(uz_mod(1:nx-1,1:nz-1))**2)
                
        open (unit=34,file='abs_u.out',status='unknown')
	    do k=1, nz-1
	    do i=1, nx-1
	    
            WRITE(34,"(1X, E18.8E3,$)")  tbplt(i,k)
        ENDDO
	    WRITE(34,*)
        ENDDO
        CLOSE(34)
                 
    endif

    enddo! end FDTD body

    END
