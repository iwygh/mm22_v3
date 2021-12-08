program midplane
implicit none
real*8, parameter::pi = 2.0d0*dasin(1.0d0)
character*1 maxustr
character*6 Nrepstr
character*6 aminstr,amaxstr,setstr,ctUstr
character*500 infile,outfile,outfilestem,dummytext,classhere
character*500 :: amin_strlist(16),amax_strlist(16)
character*500 :: amin_strlist_2021(16),amax_strlist_2021(16)
integer ct_U,line_count,max_U,dummyinteger,U_here,iu
integer Nrep,irep,Nobj,iobj,Nbin,ibin
real*8 Nsteps,amin,amax,dummyreal,ahere,xhat,yhat,zhat,time0,time1
real*8 aa,ee,ii,argperi,node,mean_anomaly,eccentric_anomaly,true_anomaly,rr,eclon,eclat
real*8 i_vm17method,Omega_vm17method,q_vm17method,p_vm17method,sum_vm17method
real*8 :: amin_list(16),amax_list(16)
real*8 :: amin_list_2021(16),amax_list_2021(16)
real*8,allocatable :: a_list(:),e_list(:),i_list(:),node_list(:),radius_list(:)
real*8,allocatable :: eclon_list(:),eclat_list(:),Omega_list(:)
real*8,allocatable :: xhat_list(:),yhat_list(:),zhat_list(:)
! from http://personal.ph.surrey.ac.uk/~phs3ps/fortweb/pre-2011/glossary/random_seed.htm
! ----- variables for portable seed setting -----
INTEGER :: i_seed
INTEGER, DIMENSION(:), ALLOCATABLE :: a_seed
INTEGER, DIMENSION(1:8) :: dt_seed
call cpu_time(time0)
! ----- end of variables for seed setting -----
! ----- Set up random seed portably -----
CALL RANDOM_SEED(size=i_seed)
ALLOCATE(a_seed(1:i_seed))
CALL RANDOM_SEED(get=a_seed)
CALL DATE_AND_TIME(values=dt_seed)
a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)
CALL RANDOM_SEED(put=a_seed)
DEALLOCATE(a_seed)
! ----- Done setting up random seed -----
Nbin = 16
Nsteps = 1000.0d0



max_U = 9
maxustr = '9'
amin_list_2021 = (/ 30.0d0,34.0d0,34.79d0,34.79d0,34.79d0,40.0d0,40.524d0,42.0d0,42.0d0,42.0d0,&
                  & 43.0d0,44.0d0,45.0d0, 45.0d0, 50.0d0, 50.0d0 /)
amax_list_2021 = (/ 150.0d0,36.0d0,40.524d0,50.0d0,150.0d0,41.0d0,42.0d0,43.0d0,45.0d0,48.0d0,&
                  & 44.0d0, 45.0d0,48.0d0,  50.0d0,80.0d0, 150.0d0 /)
amin_strlist_2021 = [ '30    ','34    ','34.79 ','34.79 ','34.79 ','40    ','40.524', &
                    & '42    ','42    ','42    ','43    ','44    ','45    ','45    ','50    ','50    ' ]
amax_strlist_2021 = [ '150   ','36    ','40.524','50    ','150   ','41    ','42    ',&
                    & '43    ','45    ','48    ','44    ','45    ','48    ','50    ','80    ','150   ' ]
infile      = 'sbdb_query_results_delcols.csv'
outfilestem = 'sbdb_query_results_delcols_'
amin_list = amin_list_2021
amax_list = amax_list_2021
amin_strlist = amin_strlist_2021
amax_strlist = amax_strlist_2021
maxustr = trim(adjustl(maxustr))
Nrepstr = trim(adjustl(Nrepstr))
infile = trim(adjustl(infile))
do ibin = 1,Nbin
    ! write(*,*) maxustr,Nrepstr
    amin = amin_list(ibin)
    amax = amax_list(ibin)
    aminstr = trim(adjustl(amin_strlist(ibin)))
    amaxstr = trim(adjustl(amax_strlist(ibin)))
    !     write(*,*) aminstr,amaxstr
    !     write(*,*) ibin,Nbin,amin,amax
    ! count non-resonant objects in semimajor axis bin
    call get_line_count(infile,line_count)
    Nobj = line_count - 1 ! header
    ct_U = 0 ! count of objects in semimajor axis bin
    open(99,file=infile)
    read(99,*) ! skip first line (header)
    do iobj = 1,Nobj
!         write(*,*) 'reading line',iobj+1,'of',Nobj+1,'in first reading round'
        read(99,*) &
            & dummyreal,    dummyreal,   dummyreal,   dummyreal,   dummyreal,&
            & dummyreal,    dummyreal,   dummyreal,   dummyreal,   dummyreal,&
            & dummyreal,    dummyreal,   dummyreal,   dummyreal,   dummyreal,&
            & dummyreal,    dummyreal,   dummyreal,   dummyreal,   dummyreal,&
            & dummyreal,    dummyreal,   dummyreal,   dummyreal,   dummyinteger,&
            & mean_anomaly, argperi,     node,        ii,          ee,&
            & dummyreal,    aa,          dummyreal,   dummyreal,   dummytext,&
            & dummytext,    classhere,   dummyreal,   dummyinteger,dummyinteger,&
            & dummyinteger, dummyinteger,rr,          xhat,        yhat,&
            & zhat,         eclat,       eclon
!             & dummytext,   dummyreal,   dummyreal,   dummyreal,   dummyreal,&
!             & dummyreal,   dummyreal,   dummyreal,   dummyreal,   dummyreal,&
!             & U_here,      dummyinteger,dummyinteger,dummyinteger,dummyinteger,&
!             & dummyinteger,dummyreal,   dummytext,   dummyreal,   mean_anomaly,&
!             & argperi,     node,        ii,          ee,          dummyreal,&
!             & aa,          dummyreal,   dummyreal,   classhere,   dummyreal,&
!             & dummyreal,   dummyreal,   dummyreal,   dummyreal,   dummyreal,&
!             & rr,          dummyreal,   eclat,       eclon,       xhat,&
!             & yhat,        zhat
        ! want to read in mean_anomaly,argperi,node,ii,ee,aa,rr,eclon,eclat,xhat,yhat,zhat
        classhere = trim(adjustl(classhere))
    !         write(*,*) trim(adjustl(classhere))
        if (((aa .ge. amin) .and. (aa .lt. amax)) .and. &
            & (trim(adjustl(classhere)) .eq.  'Classical')) then
            ct_U = ct_U + 1
    !             write(*,*) 'nonresonant object that fits the bill',ct_U,trim(adjustl(classhere))
        end if
        if (((aa .ge. amin) .and. (aa .lt. amax)) .and. &
            & (trim(adjustl(classhere)) .eq.  'Scattering')) then
            ct_U = ct_U + 1
    !             write(*,*) 'nonresonant object that fits the bill',ct_U,trim(adjustl(classhere))
        end if
        if (((aa .ge. amin) .and. (aa .lt. amax)) .and. &
            & (trim(adjustl(classhere)) .eq.  'Detached')) then
            ct_U = ct_U + 1
    !             write(*,*) 'nonresonant object that fits the bill',ct_U,trim(adjustl(classhere))
        end if
    end do
    close(99)
    !     write(*,*) 'ct_U',ct_U
    ! allocate space for barycentric orbels of satisfactory objects: a,e,i,node,radius,eclon,eclat
    if (allocated(a_list)) then
        deallocate(a_list)
    end if
    if (allocated(e_list)) then
        deallocate(e_list)
    end if
    if (allocated(i_list)) then
        deallocate(i_list)
    end if
    if (allocated(node_list)) then
        deallocate(node_list)
    end if
    if (allocated(radius_list)) then
        deallocate(radius_list)
    end if
    if (allocated(eclon_list)) then
        deallocate(eclon_list)
    end if
    if (allocated(eclat_list)) then
        deallocate(eclat_list)
    end if
    if (allocated(Omega_list)) then
        deallocate(Omega_list)
    end if
    if (allocated(xhat_list)) then
        deallocate(xhat_list)
    end if
    if (allocated(yhat_list)) then
        deallocate(yhat_list)
    end if
    if (allocated(zhat_list)) then
        deallocate(zhat_list)
    end if
    allocate( a_list(1:ct_U),e_list(1:ct_U),i_list(1:ct_U),node_list(1:ct_U),radius_list(1:ct_U) )
    allocate( eclon_list(1:ct_U),eclat_list(1:ct_U),Omega_list(1:ct_U) )
    allocate( xhat_list(1:ct_U),yhat_list(1:ct_U),zhat_list(1:ct_U) )
    ! read in file again and save a,e,i,node,radius,eclon,eclat to input lists
    iu = 0
    open(99,file=infile)
    read(99,*) ! skip first line (header)
    do iobj = 1,Nobj
!         write(*,*) 'reading line',iobj+1,'of',Nobj+1,'in second reading round'
        read(99,*) &
            & dummyreal,    dummyreal,   dummyreal,   dummyreal,   dummyreal,&
            & dummyreal,    dummyreal,   dummyreal,   dummyreal,   dummyreal,&
            & dummyreal,    dummyreal,   dummyreal,   dummyreal,   dummyreal,&
            & dummyreal,    dummyreal,   dummyreal,   dummyreal,   dummyreal,&
            & dummyreal,    dummyreal,   dummyreal,   dummyreal,   dummyinteger,&
            & mean_anomaly, argperi,     node,        ii,          ee,&
            & dummyreal,    aa,          dummyreal,   dummyreal,   dummytext,&
            & dummytext,    classhere,   dummyreal,   dummyinteger,dummyinteger,&
            & dummyinteger, dummyinteger,rr,          xhat,        yhat,&
            & zhat,         eclat,       eclon
!             & dummytext,   dummyreal,   dummyreal,   dummyreal,   dummyreal,&
!             & dummyreal,   dummyreal,   dummyreal,   dummyreal,   dummyreal,&
!             & U_here,      dummyinteger,dummyinteger,dummyinteger,dummyinteger,&
!             & dummyinteger,dummyreal,   dummytext,   dummyreal,   mean_anomaly,&
!             & argperi,     node,        ii,          ee,          dummyreal,&
!             & aa,          dummyreal,   dummyreal,   classhere,   dummyreal,&
!             & dummyreal,   dummyreal,   dummyreal,   dummyreal,   dummyreal,&
!             & rr,          dummyreal,   eclat,       eclon,       xhat,&
!             & yhat,        zhat
        ! want to read in mean_anomaly,argperi,node,ii,ee,aa,rr,eclon,eclat,xhat,yhat,zhat
        classhere = trim(adjustl(classhere))
    !         write(*,*) trim(adjustl(classhere))
        if (((aa .ge. amin) .and. (aa .lt. amax)) .and. &
            & (trim(adjustl(classhere)) .eq.  'Classical')) then
            iu = iu + 1
    !             write(*,*) 'nonresonant object that fits the bill',iu,trim(adjustl(classhere))
            mean_anomaly = pi/180.0d0 * mean_anomaly! degrees to radians
            argperi = pi/180.0d0 * argperi ! degrees to radians
            node = pi/180.0d0 * node ! degrees to radians
            ii = pi/180.0d0 * ii ! degrees to radians
            eclon = pi/180.0d0 * eclon ! degrees to radians
            eclat = pi/180.0d0 * eclat ! degrees to radians
            a_list(iu) = aa
            e_list(iu) = ee
            i_list(iu) = ii
            node_list(iu) = node
            radius_list(iu) = rr
            eclon_list(iu) = eclon
            eclat_list(iu) = eclat
        end if
        if (((aa .ge. amin) .and. (aa .lt. amax)) .and. &
            & (trim(adjustl(classhere)) .eq.  'Scattering')) then
            iu = iu + 1
    !             write(*,*) 'nonresonant object that fits the bill',iu,trim(adjustl(classhere))
            mean_anomaly = pi/180.0d0 * mean_anomaly! degrees to radians
            argperi = pi/180.0d0 * argperi ! degrees to radians
            node = pi/180.0d0 * node ! degrees to radians
            ii = pi/180.0d0 * ii ! degrees to radians
            eclon = pi/180.0d0 * eclon ! degrees to radians
            eclat = pi/180.0d0 * eclat ! degrees to radians
            a_list(iu) = aa
            e_list(iu) = ee
            i_list(iu) = ii
            node_list(iu) = node
            radius_list(iu) = rr
            eclon_list(iu) = eclon
            eclat_list(iu) = eclat
        end if
        if (((aa .ge. amin) .and. (aa .lt. amax)) .and. &
            & (trim(adjustl(classhere)) .eq.  'Detached')) then
            iu = iu + 1
    !             write(*,*) 'nonresonant object that fits the bill',iu,trim(adjustl(classhere))
            mean_anomaly = pi/180.0d0 * mean_anomaly! degrees to radians
            argperi = pi/180.0d0 * argperi ! degrees to radians
            node = pi/180.0d0 * node ! degrees to radians
            ii = pi/180.0d0 * ii ! degrees to radians
            eclon = pi/180.0d0 * eclon ! degrees to radians
            eclat = pi/180.0d0 * eclat ! degrees to radians
            a_list(iu) = aa
            e_list(iu) = ee
            i_list(iu) = ii
            node_list(iu) = node
            radius_list(iu) = rr
            eclon_list(iu) = eclon
            eclat_list(iu) = eclat
        end if
    end do
    close(99)
    !     write(*,*) 'iu',iu,'ct_U',ct_U
    ! create string for output filename
    write (ctUstr, '(I6)') ct_U
    ctUstr = trim(adjustl(ctUstr))
    outfile = trim(adjustl(outfilestem)) // &
                & '_objct' // trim(adjustl(ctUstr)) // &
                & '_amin' // trim(adjustl(aminstr)) // &
                & '_amax' // trim(adjustl(amaxstr)) // '_nominal' // '.txt'
    outfile = trim(adjustl(outfile))
    write(*,*) outfile
    open(99,file=outfile)
    write(99,*) 'i_vm ','Omega_vm ','q_vm ','p_vm ','sum_vm '
    ! compute the mean plane of the nominal objects in the bins
    call nominal_mean_plane(amin,amax,Nsteps,ct_U,a_list,e_list,i_list,node_list,&
                & eclon_list,eclat_list,radius_list,irep,Nrep,&
                & xhat_list,yhat_list,zhat_list,& !inputs
                & i_vm17method,Omega_vm17method,q_vm17method,p_vm17method,sum_vm17method) ! outputs
    ! write the results to output file
    write(99,*) i_vm17method,Omega_vm17method,q_vm17method,p_vm17method,sum_vm17method
    close(99)
end do
stop
end program







subroutine get_line_count(infile,& ! inputs
                  & line_count ) ! outputs
implicit none
integer line_count,io
character(len=*) infile
open(99,file=infile,iostat=io)
if (io/=0) stop 'get_line_count: Cannot open file!'
line_count = 0
do
    read(99,*,iostat=io)
    if (io/=0) exit
    line_count = line_count + 1
end do
close(99)
return
end subroutine

subroutine fit_plane_qp(Nobj,Nsteps,vx_list,vy_list,vz_list,q_list,p_list,i_list,Omega_list, & ! inputs
                    & i_qp,Omega_qp,q_qp,p_qp,sum_qp) ! outputs
implicit none
integer Nobj,i
real*8 Nsteps,i_qp,Omega_qp,q_qp,p_qp,sum_qp,q_min,q_max,dq,p_min,p_max,dp, &
        & qnom,pnom,rq,rp,nx,ny,nz,q_here,p_here,sum_here,term,sin_i_qp
real*8 vx_list(Nobj),vy_list(Nobj),vz_list(Nobj),i_list(Nobj),Omega_list(Nobj), &
                    & q_list(Nobj),p_list(Nobj)
real*8, parameter::pi = 2.0d0*dasin(1.0d0)
q_min = minval(q_list)
q_max = maxval(q_list)
dq = (q_max-q_min)/Nsteps
p_min = minval(p_list)
p_max = maxval(p_list)
dp = (p_max-p_min)/Nsteps
sum_qp = 1.0d9 ! a large number
qnom = q_min
pnom = p_min
do while (qnom < q_max)
    q_here = qnom
    do while (pnom < p_max)
        p_here = pnom
        nx = p_here
        ny = -q_here
        nz = dsqrt(1.0d0-nx*nx-ny*ny) ! assume i<90deg for midplane
        sum_here = 0.0d0
        do i = 1,Nobj
            term = vx_list(i)*nx + vy_list(i)*ny + vz_list(i)*nz
            sum_here = sum_here + abs(term)
        end do
        if (sum_here .le. sum_qp) then
            q_qp = q_here
            p_qp = p_here
            sum_qp = sum_here
            sin_i_qp = dsqrt(q_qp*q_qp+p_qp*p_qp)
            i_qp = dasin(sin_i_qp) ! assume i<90deg for midplane
            Omega_qp = datan2(p_qp/sin_i_qp,q_qp/sin_i_qp)
            Omega_qp = Omega_qp*180.0d0/pi
            if(Omega_qp<0.0d0) Omega_qp=Omega_qp+360.0d0
            i_qp = i_qp*180.0d0/pi
        end if
        pnom = pnom + dp
    end do
    qnom = qnom + dq
    pnom = p_min
end do
sum_qp = sum_qp / Nobj
return
end subroutine

subroutine vm17method2(amin,amax,Nsteps,Nobj,a_list,e_list,i_list,Omega_list,&
                    & eclon_list,eclat_list,rbaryau_list,irep,Nrep,&
                    & xhat_list,yhat_list,zhat_list,& !inputs
                    & i_vm17,Omega_vm17,q_vm17,p_vm17,sum_vm17) ! outputs
implicit none
integer Nobj,i,cnt,match,j,maxiter,irep,Nrep
real*8, parameter::pi = 2.0d0*dasin(1.0d0)
real*8 i_vm17,Omega_vm17,q_vm17,p_vm17,dummy,hx,hy,hz,xhat,yhat,zhat, &
    & i_qp,Omega_qp,q_qp,p_qp,sum_qp,sigma1,sigma2,inc,var_q,var_p, &
    & std_q,std_p,std,mean_q,mean_p,mean,sum_q,sum_p,w,temp,temp2,x1,x2, &
    & q,p,node,aperi,ma,at,amin,amax,et,GM,x,y,z,vx,vy,vz,r,rp,eclat,eclon, &
    & rtarget,diff_lat,diff_lon,dr,q0,p0,Nsteps,sum_vm17
real*8 i_list(Nobj),Omega_list(Nobj),eclon_list(Nobj),eclat_list(Nobj), &
    & vx_list(Nobj),vy_list(Nobj),vz_list(Nobj),a_list(Nobj),e_list(Nobj), &
    & rbaryau_list(Nobj),vx_list_2(Nobj),vy_list_2(Nobj),vz_list_2(Nobj), &
    & i_list_2(Nobj),Omega_list_2(Nobj),q_list(Nobj),p_list(Nobj),q_list_2(Nobj),p_list_2(Nobj),&
    & xhat_list(Nobj),yhat_list(Nobj),zhat_list(Nobj)
GM = 1.0d0
maxiter = 10000
! write(*,*) "starting vm17method"
do i = 1,Nobj
    hx = dsin(i_list(i))*dsin(Omega_list(i))
    hy = -dsin(i_list(i))*dcos(Omega_list(i))
    hz = dcos(i_list(i))
    xhat = dcos(eclat_list(i))*dcos(eclon_list(i))
    yhat = dcos(eclat_list(i))*dsin(eclon_list(i))
    zhat = dsin(eclat_list(i))
    vx_list(i) = hy*zhat - hz*yhat
    vy_list(i) = hz*xhat - hx*zhat
    vz_list(i) = hx*yhat - hy*xhat
!     vx_list(i) = hy*zhat_list(i) - hz*yhat_list(i)
!     vy_list(i) = hz*xhat_list(i) - hx*zhat_list(i)
!     vz_list(i) = hx*yhat_list(i) - hy*xhat_list(i)
    q_list(i) = -hy
    p_list(i) = hx
end do
call fit_plane_qp(Nobj,Nsteps,vx_list,vy_list,vz_list,q_list,p_list,i_list,Omega_list, & ! inputs
                & i_qp,Omega_qp,q_qp,p_qp,sum_qp) ! outputs
q0 = q_qp
p0 = p_qp
! get estimate of Gaussian qp stdev of input population
! assuming q,p are independent and identically distributed,
! ie q and p scatter in a circle with zero eccentricity or correlation
sum_q = 0.0d0
sum_p = 0.0d0
do i = 1,Nobj
    sum_q = sum_q + q_list(i)
    sum_p = sum_p + p_list(i)
end do
mean_q = sum_q / Nobj
mean_p = sum_p / Nobj
var_q = 0.0d0
var_p = 0.0d0
do i = 1,Nobj
    var_q = var_q + (q_list(i)-mean_q)**2
    var_p = var_p + (p_list(i)-mean_p)**2
end do
var_q = var_q / (Nobj-1)
var_p = var_p / (Nobj-1)
std_q = dsqrt(var_q)
std_p = dsqrt(var_p)
do i = 1,Nobj
    cnt = 0
    match = 0
!         write(*,*) "Trying to generate synthetic object no. ",i
!         call sleep(2)
    do while ((match .eq. 0) .and. (cnt .le. maxiter))
        cnt = cnt + 1
        if (cnt .eq. (maxiter-1) ) then
            write(*,*) "maxiter",maxiter," reached for object",i,"of",Nobj,&
                & "on rep",irep,"of",Nrep,"for amin =",amin,"amax =",amax
        end if
!            write(*,*) "Try",cnt,"for synthetic object",i,"of",NObj,&
!                & "on rep",irep,"of",Nrep,"for amin =",amin,"amax =",amax
        ! pulling inclination components from the Gaussian
        inc = 1.5d0
        do while (inc .ge. 1d0)
            w = 1.5d0
            do while (w .ge. 1d0)
                call random_number(temp)
                call random_number(temp2)
                x1 = 2.0d0*temp - 1.0d0
                x2 = 2.0d0*temp2 - 1.0d0
                w = x1*x1 + x2*x2
            end do
            w = dsqrt( (-2.0d0*dlog(w)) / w )
            temp = x1*w
            temp2 = x2*w
            ! assigning inclination components relative to the forced plane
            q = q0 + temp*std_q
            p = p0 + temp2*std_p
            inc = dsqrt(q*q+p*p)
        end do
        inc = dasin(inc)
        node = datan2(p,q)
        if (node .lt. 0.0d0) node = node + 2.0d0*pi
        ! randomly assigning argument of perihelion
        call random_number(temp)
        aperi = temp*2.0d0*pi
        ! randomly assigning mean anomaly
        call random_number(temp)
        ma = temp*2.0d0*pi
!             randomly assigning semimajor axis from amin to amax
!             call random_number(temp)
!             at = amin + temp*(amax-amin)
        ! randomly select one of the observed objects
        j = 0
        do while (j .eq. 0 .or. j .gt. Nobj)
            call random_number(temp)
            j = nint(temp*Nobj+0.5d0)
        end do
        ! then use that object's eccentricity and semimajor axis (fuzzed a little)
        call random_number(temp)
        et = e_list(j) + (0.5d0-temp)*0.05d0*e_list(j)
        call random_number(temp)
        at = a_list(j) + (0.5d0-temp)*0.01d0*a_list(j)
        ! now convert that fully specified orbit to cartesian coordinates
        call el2xv(GM,at,et,inc,node,aperi,ma, & ! inputs
            & x,y,z,vx,vy,vz) ! outputs
        r = x*x + y*y + z*z
        r = dsqrt(r)
        rp = x*x + y*y
        rp = dsqrt(rp)
        ! then convert to sky position
!             eclat = datan2(z,rp)
        eclat = dasin(z/r)
        eclon = datan2(y,x)
        if (eclon .lt. 0.0d0) eclon = eclon + 2.0d0*pi
        ! and check to see if it's near the observed object
        rtarget = rbaryau_list(j)
        dr = dabs(r-rtarget)/r
        diff_lat = dabs(eclat*pi/180.0d0 - eclat_list(i)*pi/180.0d0) ! easier to compare in degrees than radians
        diff_lon = dabs(eclon*pi/180.0d0 - eclon_list(i)*pi/180.0d0) ! easier to compare in degrees than radians
        temp = dabs(diff_lon-360.0d0)
        if (temp .lt. diff_lon) diff_lon = temp ! a difference of 359 degrees is the same as a difference of 1 degree
        ! match if eclon diff 5 deg, eclat diff 1 deg, distance diff 10%
!             if ( (diff_lon .lt. 5.0d0) .and. (diff_lat .le. 1.0d0) .and. (dr .le. 0.1d0) ) match = 1
        if ( (diff_lon .le. 5.0d0) .and. (diff_lat .le. 1.0d0)) match = 1
    end do
    ! we found a match, compute vx,vy,vz,i,Omega
    hx = dsin(inc)*dsin(node)
    hy = -dsin(inc)*dcos(node)
    hz = dcos(inc)
    xhat = dcos(eclat)*dcos(eclon)
    yhat = dcos(eclat)*dsin(eclon)
    zhat = dsin(eclat)
    vx_list_2(i) = hy*zhat - hz*yhat
    vy_list_2(i) = hz*xhat - hx*zhat
    vz_list_2(i) = hx*yhat - hy*xhat
    i_list_2(i) = inc
    Omega_list_2(i) = node
    q_list_2(i) = dsin(inc)*dcos(node)
    p_list_2(i) = dsin(inc)*dsin(node)
end do
! get midplane estimate of synthetic population
! write(*,*) "computing midplane for synthetic population",irep,"of",Nrep
call fit_plane_qp(Nobj,Nsteps,vx_list_2,vy_list_2,vz_list_2,q_list_2,p_list_2,i_list_2,Omega_list_2, & ! inputs
                & i_vm17,Omega_vm17,q_vm17,p_vm17,sum_vm17) ! outputs
return
end subroutine


subroutine nominal_mean_plane(amin,amax,Nsteps,Nobj,a_list,e_list,i_list,Omega_list,&
                    & eclon_list,eclat_list,rbaryau_list,irep,Nrep,&
                    & xhat_list,yhat_list,zhat_list,& !inputs
                    & i_vm17,Omega_vm17,q_vm17,p_vm17,sum_vm17) ! outputs
! modified from vm17method2 to do only the nominal mean plane without all the random stuff
implicit none
integer Nobj,i,cnt,match,j,maxiter,irep,Nrep
real*8, parameter::pi = 2.0d0*dasin(1.0d0)
real*8 i_vm17,Omega_vm17,q_vm17,p_vm17,dummy,hx,hy,hz,xhat,yhat,zhat, &
    & i_qp,Omega_qp,q_qp,p_qp,sum_qp,sigma1,sigma2,inc,var_q,var_p, &
    & std_q,std_p,std,mean_q,mean_p,mean,sum_q,sum_p,w,temp,temp2,x1,x2, &
    & q,p,node,aperi,ma,at,amin,amax,et,GM,x,y,z,vx,vy,vz,r,rp,eclat,eclon, &
    & rtarget,diff_lat,diff_lon,dr,q0,p0,Nsteps,sum_vm17
real*8 i_list(Nobj),Omega_list(Nobj),eclon_list(Nobj),eclat_list(Nobj), &
    & vx_list(Nobj),vy_list(Nobj),vz_list(Nobj),a_list(Nobj),e_list(Nobj), &
    & rbaryau_list(Nobj),vx_list_2(Nobj),vy_list_2(Nobj),vz_list_2(Nobj), &
    & i_list_2(Nobj),Omega_list_2(Nobj),q_list(Nobj),p_list(Nobj),q_list_2(Nobj),p_list_2(Nobj),&
    & xhat_list(Nobj),yhat_list(Nobj),zhat_list(Nobj)
GM = 1.0d0
maxiter = 10000
! write(*,*) "starting vm17method"
do i = 1,Nobj
    hx = dsin(i_list(i))*dsin(Omega_list(i))
    hy = -dsin(i_list(i))*dcos(Omega_list(i))
    hz = dcos(i_list(i))
    xhat = dcos(eclat_list(i))*dcos(eclon_list(i))
    yhat = dcos(eclat_list(i))*dsin(eclon_list(i))
    zhat = dsin(eclat_list(i))
    vx_list(i) = hy*zhat - hz*yhat
    vy_list(i) = hz*xhat - hx*zhat
    vz_list(i) = hx*yhat - hy*xhat
!     vx_list(i) = hy*zhat_list(i) - hz*yhat_list(i)
!     vy_list(i) = hz*xhat_list(i) - hx*zhat_list(i)
!     vz_list(i) = hx*yhat_list(i) - hy*xhat_list(i)
    q_list(i) = -hy
    p_list(i) = hx
end do
call fit_plane_qp(Nobj,Nsteps,vx_list,vy_list,vz_list,q_list,p_list,i_list,Omega_list, & ! inputs
                & i_qp,Omega_qp,q_qp,p_qp,sum_qp) ! outputs
i_vm17 = i_qp
Omega_vm17 = Omega_qp
q_vm17 = q_qp
p_vm17 = p_qp
return
end subroutine







subroutine vm17method(amin,amax,irep,Nrep,infile,Nsteps,Nobj,& ! inputs
                & i_vm17,Omega_vm17,q_vm17,p_vm17,sum_qp) ! outputs
implicit none
character(len=*) infile
integer Nobj,i,cnt,match,j,irep,Nrep,maxiter
real*8, parameter::pi = 2.0d0*dasin(1.0d0)
real*8 i_vm17,Omega_vm17,q_vm17,p_vm17,dummy,hx,hy,hz,xhat,yhat,zhat, &
    & i_qp,Omega_qp,q_qp,p_qp,sum_qp,sigma1,sigma2,inc,var_q,var_p, &
    & std_q,std_p,std,mean_q,mean_p,mean,sum_q,sum_p,w,temp,temp2,x1,x2, &
    & q,p,node,aperi,ma,at,amin,amax,et,GM,x,y,z,vx,vy,vz,r,rp,eclat,eclon, &
    & rtarget,diff_lat,diff_lon,dr,q0,p0,Nsteps
real*8 i_list(Nobj),Omega_list(Nobj),eclon_list(Nobj),eclat_list(Nobj), &
    & vx_list(Nobj),vy_list(Nobj),vz_list(Nobj),a_list(Nobj),e_list(Nobj), &
    & rbaryau_list(Nobj),vx_list_2(Nobj),vy_list_2(Nobj),vz_list_2(Nobj), &
    & i_list_2(Nobj),Omega_list_2(Nobj),q_list(Nobj),p_list(Nobj),q_list_2(Nobj),p_list_2(Nobj)
GM = 1.0d0
maxiter = 10000
write(*,*) "starting vm17method"
! if (Nobj>0) then
    open(100,file=infile)
    read(100,*) ! skip first line (header)
    do i = 1,Nobj ! read in data
        read(100,*) dummy,a_list(i),e_list(i),dummy,dummy,dummy,dummy, &
                & i_list(i),dummy,Omega_list(i),eclon_list(i),eclat_list(i), &
                & rbaryau_list(i)
        if (eclon_list(i)<0.0d0) eclon_list(i) = eclon_list(i) + 2.0d0*pi
        hx = dsin(i_list(i))*dsin(Omega_list(i))
        hy = -dsin(i_list(i))*dcos(Omega_list(i))
        hz = dcos(i_list(i))
        xhat = dcos(eclat_list(i))*dcos(eclon_list(i))
        yhat = dcos(eclat_list(i))*dsin(eclon_list(i))
        zhat = dsin(eclat_list(i))
        vx_list(i) = hy*zhat - hz*yhat
        vy_list(i) = hz*xhat - hx*zhat
        vz_list(i) = hx*yhat - hy*xhat
        q_list(i) = -hy
        p_list(i) = hx
    end do
    close(100)
    ! get estimate of semimajor axis bounds of input population
!     amin = minval(a_list)
!     amax = maxval(a_list)
    ! get midplane estimate of input population
    call fit_plane_qp(Nobj,Nsteps,vx_list,vy_list,vz_list,q_list,p_list,i_list,Omega_list, & ! inputs
                    & i_qp,Omega_qp,q_qp,p_qp,sum_qp) ! outputs
!     q0 = dsin(i_qp*pi/180.0d0)*dcos(Omega_qp*pi/180.0d0)
!     p0 = dsin(i_qp*pi/180.0d0)*dsin(Omega_qp*pi/180.0d0)
    q0 = q_qp
    p0 = p_qp
    ! get estimate of Gaussian qp stdev of input population
    ! assuming q,p are independent and identically distributed,
    ! ie q and p scatter in a circle with zero eccentricity or correlation
    sum_q = 0.0d0
    sum_p = 0.0d0
    do i = 1,Nobj
        sum_q = sum_q + q_list(i)
        sum_p = sum_p + p_list(i)
    end do
    mean_q = sum_q / Nobj
    mean_p = sum_p / Nobj
    var_q = 0.0d0
    var_p = 0.0d0
    do i = 1,Nobj
        var_q = var_q + (q_list(i)-mean_q)**2
        var_p = var_p + (p_list(i)-mean_p)**2
    end do
    var_q = var_q / (Nobj-1)
    var_p = var_p / (Nobj-1)
!     std = 0.5d0 * (var_q + var_p)
!     std = dsqrt(std)
    std_q = dsqrt(var_q)
    std_p = dsqrt(var_p)
!     std = (std_q+std_p)*0.5d0
    ! draw synthetic population
    do i = 1,Nobj
        cnt = 0
        match = 0
!         write(*,*) "Trying to generate synthetic object no. ",i
!         call sleep(2)
        do while ((match .eq. 0) .and. (cnt .le. maxiter))
            cnt = cnt + 1
            if (cnt .eq. (maxiter-1) ) then
                write(*,*) "maxiter",maxiter," reached for object",i,"of",Nobj,&
                    & "on rep",irep,"of",Nrep,"for amin =",amin,"amax =",amax
            end if
!            write(*,*) "Try",cnt,"for synthetic object",i,"of",NObj,&
!                & "on rep",irep,"of",Nrep,"for amin =",amin,"amax =",amax
            ! pulling inclination components from the Gaussian
            inc = 1.5d0
            do while (inc .ge. 1d0)
                w = 1.5d0
                do while (w .ge. 1d0)
                    call random_number(temp)
                    call random_number(temp2)
                    x1 = 2.0d0*temp - 1.0d0
                    x2 = 2.0d0*temp2 - 1.0d0
                    w = x1*x1 + x2*x2
                end do
                w = dsqrt( (-2.0d0*dlog(w)) / w )
                temp = x1*w
                temp2 = x2*w
                ! assigning inclination components relative to the forced plane
                q = q0 + temp*std_q
                p = p0 + temp2*std_p
                inc = dsqrt(q*q+p*p)
            end do
            inc = dasin(inc)
            node = datan2(p,q)
            if (node .lt. 0.0d0) node = node + 2.0d0*pi
            ! randomly assigning argument of perihelion
            call random_number(temp)
            aperi = temp*2.0d0*pi
            ! randomly assigning mean anomaly
            call random_number(temp)
            ma = temp*2.0d0*pi
!             randomly assigning semimajor axis from amin to amax
!             call random_number(temp)
!             at = amin + temp*(amax-amin)
            ! randomly select one of the observed objects
            j = 0
            do while (j .eq. 0 .or. j .gt. Nobj)
                call random_number(temp)
                j = nint(temp*Nobj+0.5d0)
            end do
            ! then use that object's eccentricity and semimajor axis (fuzzed a little)
            call random_number(temp)
            et = e_list(j) + (0.5d0-temp)*0.05d0*e_list(j)
            call random_number(temp)
            at = a_list(j) + (0.5d0-temp)*0.01d0*a_list(j)
            ! now convert that fully specified orbit to cartesian coordinates
            call el2xv(GM,at,et,inc,node,aperi,ma, & ! inputs
                & x,y,z,vx,vy,vz) ! outputs
            r = x*x + y*y + z*z
            r = dsqrt(r)
            rp = x*x + y*y
            rp = dsqrt(rp)
            ! then convert to sky position
!             eclat = datan2(z,rp)
            eclat = dasin(z/r)
            eclon = datan2(y,x)
            if (eclon .lt. 0.0d0) eclon = eclon + 2.0d0*pi
            ! and check to see if it's near the observed object
            rtarget = rbaryau_list(j)
            dr = dabs(r-rtarget)/r
            diff_lat = dabs(eclat*pi/180.0d0 - eclat_list(i)*pi/180.0d0) ! easier to compare in degrees than radians
            diff_lon = dabs(eclon*pi/180.0d0 - eclon_list(i)*pi/180.0d0) ! easier to compare in degrees than radians
            temp = dabs(diff_lon-360.0d0)
            if (temp .lt. diff_lon) diff_lon = temp ! a difference of 359 degrees is the same as a difference of 1 degree
            ! match if eclon diff 5 deg, eclat diff 1 deg, distance diff 10%
!             if ( (diff_lon .lt. 5.0d0) .and. (diff_lat .le. 1.0d0) .and. (dr .le. 0.1d0) ) match = 1
            if ( (diff_lon .le. 5.0d0) .and. (diff_lat .le. 1.0d0)) match = 1
        end do
        ! we found a match, compute vx,vy,vz,i,Omega
        hx = dsin(inc)*dsin(node)
        hy = -dsin(inc)*dcos(node)
        hz = dcos(inc)
        xhat = dcos(eclat)*dcos(eclon)
        yhat = dcos(eclat)*dsin(eclon)
        zhat = dsin(eclat)
        vx_list_2(i) = hy*zhat - hz*yhat
        vy_list_2(i) = hz*xhat - hx*zhat
        vz_list_2(i) = hx*yhat - hy*xhat
        i_list_2(i) = inc
        Omega_list_2(i) = node
        q_list_2(i) = dsin(inc)*dcos(node)
        p_list_2(i) = dsin(inc)*dsin(node)
    end do
    ! get midplane estimate of synthetic population
    ! write(*,*) "computing midplane for synthetic population",irep,"of",Nrep
    call fit_plane_qp(Nobj,Nsteps,vx_list_2,vy_list_2,vz_list_2,q_list_2,p_list_2,i_list_2,Omega_list_2, & ! inputs
                    & i_vm17,Omega_vm17,q_vm17,p_vm17,sum_qp) ! outputs
! else
!     i_vm17 = -99999
!     Omega_vm17 = -99999
!     q_vm17 = -99999
!     p_vm17 = -99999
!     sum_qp = -99999
! end if
return
end subroutine


!******************************************************************************
!                             SUBROUTINE EL2XV
!******************************************************************************
! Computes cartesian positions and velocities given central mass, orbital
! elements, and ialpha
! Arguments:
!   GM: G times central mass
!   \\\\\ taken out!  now only good for ellipses
!   \\\\\\  ialpha: integer for conic section type ( =+1 for hyperbola, 0 for parabola,
!   \\\\\\\        and -1 for ellipse)
!   a : semi-major axis (or pericentric distance if a parabola),
!   e : eccentricity
!   inc : inclination
!   capom : longitude of ascending node , often called capital OMEGA
!   omega : argument of perihelion (from ascending node)
!   capmnq : either M, N or Q.
! ALGORITHM:  See Fitzpatrick "Principles of Cel. Mech."
! AUTHOR:  M. Duncan.
! DATE WRITTEN:  May 11, 1992.
! REVISIONS: May 26 - now use better Kepler solver for ellipses
!   and hyperbolae called EHYBRID.F and FHYBRID.F
! Last change: 19 Feb 1998 // altered Dec 14,2009
!******************************************************************************
subroutine el2xv(GM,a,e,inc,capom,omega,capmnq,x,y,z,vx,vy,vz)
implicit none
!..............................................................................
! arguments
real*8 GM,a,e,inc,capom,omega,capmnq
real*8 x,y,z,vx,vy,vz
!..............................................................................
! internal variables
real*8 ehybrid,cape,fhybrid,capf,zget,zpara
real*8 sp,cp,so,co,si,ci
real*8 d11,d12,d13,d21,d22,d23
real*8 scap,ccap,shcap,chcap
real*8 sqe,sqgma,xfac1,xfac2,ri,vfac1,vfac2
!------------------------------------------------------------------------------
! Generate rotation matrices (on p. 42 of Fitzpatrick)
call scget(omega,sp,cp)
call scget(capom,so,co)
call scget(inc,si,ci)
d11 = cp*co - sp*so*ci
d12 = cp*so + sp*co*ci
d13 = sp*si
d21 = -sp*co - cp*so*ci
d22 = -sp*so + cp*co*ci
d23 = cp*si
! Get the other quantities depending on orbit type ( i.e. IALPHA)
! only valid for ellipses
cape = ehybrid(e,capmnq)
call scget(cape,scap,ccap)
sqe = dsqrt(1.d0 -e*e)
sqgma = dsqrt(GM*a)
xfac1 = a*(ccap - e)
xfac2 = a*sqe*scap
ri = 1.d0/(a*(1.d0 - e*ccap))
vfac1 = -ri * sqgma * scap
vfac2 = ri * sqgma * sqe * ccap
x =  d11*xfac1 + d21*xfac2
y =  d12*xfac1 + d22*xfac2
z =  d13*xfac1 + d23*xfac2
vx = d11*vfac1 + d21*vfac2
vy = d12*vfac1 + d22*vfac2
vz = d13*vfac1 + d23*vfac2
return
end subroutine

!******************************************************************************
real*8 FUNCTION EHYBRID(E,M)
!******************************************************************************
!       PURPOSE:  Solves Kepler's eqn.
!       ARGUMENTS:  Input is eccentricity E and mean anomaly M.
!       RETURNS the eccentric anomaly EHYBRID
!       ALGORITHM: For e < 0.18 uses fast routine ESOLMD
!	         For larger e but less than 0.8, uses EGET
!	         For e > 0.8 uses EHIE
!       REMARKS: Only EHIE brings M and E into range (0,TWOPI)
!       AUTHOR: M. Duncan
!       DATE WRITTEN: May 25,1992.
!******************************************************************************
implicit none
real*8 e,m,esolmd,eget,ehie
if(e .lt. 0.18d0) then
  EHYBRID = esolmd(e,m)
else
  if( e .le. 0.8d0) then
     EHYBRID = eget(e,m)
  else
     EHYBRID = ehie(e,m)
  endif
endif
return
end function

!******************************************************************************
real*8 FUNCTION ESOLMD(E,M)
!******************************************************************************
!       PURPOSE:  Solves Kepler's eqn.
!       ARGUMENTS:  Input is eccentricity E and mean anomaly EM.
!       RETURNs the eccentric anomaly ESOLVE
!       ALGORITHM: Some sort of quartic convergence from Wisdom.
!       REMARKS: Only good for small eccentricity since it only
!         iterates once. (good for planet orbits)
!      	  also does not put M between 0. and 2*pi
!       INCLUDES: needs SCGET.F
!       AUTHOR: M. Duncan
!       DATE WRITTEN: May 7, 1992.
!******************************************************************************
implicit none
real*8 x,e,m,sm,cm,sx,cx
real*8 es,ec,f,fp,fpp,fppp,dx
call scget(m,sm,cm)
x = m + e*sm*( 1.d0 + e*( cm + e*( 1.d0 -1.5d0*sm*sm)))
call scget(x,sx,cx)
es = e*sx
ec = e*cx
f = x - es  - m
fp = 1.d0 - ec
fpp = es
fppp = ec
dx = -f/fp
dx = -f/(fp + dx*fpp/2.d0)
dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0)
esolmd = x + dx
return
end function

!******************************************************************************
real*8 FUNCTION EGET(E,M)
!******************************************************************************
!       PURPOSE:  Solves Kepler's eqn.
!       ARGUMENTS:  Input is eccentricity E and mean anomaly M.
!       RETURNs the eccentric anomaly EGET
!       ALGORITHM: Quartic convergence from Danby
!       REMARKS: For results very near roundoff, give it M between
!           0 and 2*pi. One can condition M before calling EGET
!           by calling my double precision function MOD2PI(M).
!           This is not done within the routine to speed it up
!           and because it works fine even for large M.
!       AUTHOR: M. Duncan
!       DATE WRITTEN: May 7, 1992.
!       REVISIONS: May 21, 1992.  Now have it go through EXACTLY two
!           iterations with the premise that it will only be called if
!           we have an ellipse with e between 0.15 and 0.8
!******************************************************************************
! MAY 21 : FOR e < 0.18 use ESOLMD for speed and sufficient accuracy
! MAY 21 : FOR e > 0.8 use EHIE - this one may not converge fast enough.
implicit none
real*8 x,e,m,sm,cm,sx,cx
real*8 es,ec,f,fp,fpp,fppp,dx
call scget(m,sm,cm)
! begin with a guess accurate to order ecc**3
x = m + e*sm*( 1.d0 + e*( cm + e*( 1.d0 -1.5d0*sm*sm)))
! go through one iteration for improved estimate
call scget(x,sx,cx)
es = e*sx
ec = e*cx
f = x - es  - m
fp = 1.d0 - ec
fpp = es
fppp = ec
dx = -f/fp
dx = -f/(fp + dx*fpp/2.d0)
dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0)
eget = x + dx
! Do another iteration.
! For m between 0 and 2*pi this seems to be enough to get near roundoff
! error for eccentricities between 0 and 0.8
x = eget
call scget(x,sx,cx)
es = e*sx
ec = e*cx
f = x - es  - m
fp = 1.d0 - ec
fpp = es
fppp = ec
dx = -f/fp
dx = -f/(fp + dx*fpp/2.d0)
dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0)
eget = x + dx
return
end function

!******************************************************************************
real*8 FUNCTION EHIE(E,M)
!******************************************************************************
!       PURPOSE:  Solves Kepler's eqn.
!       ARGUMENTS:  Input is eccentricity E and mean anomaly M.
!       RETURNs the eccentric anomaly EHIE
!       ALGORITHM: Use Danby's quartic for 3 iterations.
!                Eqn. is f(x) = x - e*sin(x+M).
!                Note that E = x + M first guess is very good for e near 1.
!	         Need to first get M between 0. and PI and use symmetry
!                to return right answer if M between PI and 2PI
!       REMARKS: Modifies M so that both E and M are in range (0,TWOPI)
!       AUTHOR: M. Duncan
!       DATE WRITTEN: May 25,1992.
!******************************************************************************
implicit none
integer iflag,nper,niter
real*8 e,m
real*8 dx,x,sa,ca,esa,eca,f,fp
real*8, parameter:: PI = 3.1415926535897932d0
real*8, parameter::TWOPI = 2.d0*PI
real*8, parameter ::TOL = 3.d-15
integer, parameter::NPLMAX = 3
!	Bring M into the range (0,TWOPI), and if the result is greater than PI,
!	solve for (TWOPI - M).
iflag = 0
nper = m/TWOPI
m = m - nper*TWOPI
if (m .lt. 0.d0) m = m + TWOPI
if (m.gt.PI) then
   m = TWOPI - m
   iflag = 1
endif
!	Make a first guess that works well for e near 1.
x = (6.d0*m)**(1.d0/3.d0) - m
niter =0
!	Iteration loop
do niter =1,NPLMAX
    call scget(x + m,sa,ca)
    esa = e*sa
    eca = e*ca
    f = x - esa
    fp = 1.d0 -eca
    dx = -f/fp
    dx = -f/(fp + 0.5d0*dx*esa)
    dx = -f/(fp + 0.5d0*dx*(esa+0.3333333333333333d0*eca*dx))
!	    write(6,*) niter,m,dx
    x = x + dx
enddo
EHIE = m + x
if (iflag.eq.1) then
  EHIE = TWOPI - EHIE
  m = TWOPI - m
endif
return
end function

!******************************************************************************
!	                        SUBROUTINE SCGET
!******************************************************************************
!       PURPOSE:  Given an angle (in RADIANS), efficiently compute sin and cos.
!       ARGUMENTS:  Input is a real*8 ANGLE in radians.
!                 Returned are SX=sin(ANGLE) and CX=cos(ANGLE)
!       REMARKS: The HP 700 series won't return correct answers for sin
!         and cos if the angle is bigger than 3e7. We first reduce it
!         to the range [0,2pi) and use the sqrt rather than cos (it's faster)
!         BE SURE THE ANGLE IS IN RADIANS - NOT DEGREES!
!       AUTHOR:  M. Duncan.
!       DATE WRITTEN:  May 6, 1992.
!       Last change: 19 Feb 1998 (RM)
!******************************************************************************
subroutine scget(angle,sx,cx)
implicit none
integer nper
real*8 angle,sx,cx
real*8 x
real*8, parameter:: PI = 3.1415926535897932d0
real*8, parameter:: TWOPI = 2.d0*PI
real*8, parameter:: PIBY2 = 0.5d0*PI
real*8, parameter:: PI3BY2 = 1.5d0*PI
nper = angle/TWOPI
x = angle - nper*TWOPI
if(x.lt.0.d0) x = x + TWOPI
sx = dsin(x)
cx= dsqrt(1.d0 - sx*sx)
if( x .gt. PIBY2 .and. x .lt.PI3BY2) cx = -cx
return
end subroutine

!******************************************************************************
! 			       SUBROUTINE XV2EL
!******************************************************************************
!     PURPOSE:  Given the cartesian position and velocity of an orbit,
!       compute the osculating orbital elements.
!     INPUTS: position (X,Y,Z), velocity (VX,VY,VZ), and central GM
!     RETURNED:
!       ialpha: an integer: (-1 for ellipse, 0 for parabola, +1 for hyperbola)
!       a:  semi-major axis for ellipse, pericentric distance for parabola,
!           (-2*GM/energy) for hyperbola.
!       e: eccentricity
!       inc: inclination
!       node: longitude of the ascending node
!       peri: argument of perihelion
!	capmnq: mean anomaly M for an ellipse, Q for a parabola or
!               N for a hyperbola (in Fitzpatrick's notation).
!       All angles are in radians.
!     ALGORITHM: See e.g. p.70 of Fitzpatrick's "Priciples of Cel. Mech."
!     REMARKS:
!     1. If the inclination INC is less than TINY, we arbitrarily choose the
!        longitude of the ascending node NODE to be 0.0 (so the ascending
!        node is then along the X axis).
!     2. If the eccentricity, E, is less than SQRT(TINY), we arbitrarily
!        choose the argument of perihelion PERI to be 0.
!     AUTHOR: M. Duncan; May 8,1992.
! Last change: 28 Jan 1998 (RM)
!******************************************************************************
subroutine xv2el(GM,x,y,z,vx,vy,vz,a,e,inc,node,peri,capmnq)
implicit none
integer ialpha
real*8 x,y,z,vx,vy,vz
real*8 GM,a,e,inc,node,peri,capmnq
real*8 hx,hy,hz,h2,h,r,v,v2,vdotr,energy
real*8 es,ec,cw,sw,w,u
real*8 fac,cape,capf,tmpf
real*8, PARAMETER:: tiny = 2.D-15
real*8, parameter:: pi=3.1415926535897932d0,tpi=2.d0*pi
! Compute the angular momentum H, and thereby the inclination INC.
hx = y*vz - z*vy
hy = z*vx - x*vz
hz = x*vy - y*vx
h2 = hx*hx + hy*hy +hz*hz
h  = dsqrt(h2)
inc = dacos(hz/h)
! Compute longitude of ascending node CAPOM and the argument of latitude u
fac = dsqrt(hx*hx + hy*hy)/h
if(fac.lt. TINY ) then
  node = 0.d0
  u = datan2(y,x)
  if(dabs(inc - pi).lt. 10.d0*TINY) u = -u
else
  node = datan2(hx,-hy)
  u = datan2(z/dsin(inc), x*dcos(node) + y*dsin(node))
endif
if(node .lt. 0.d0) node = node + tpi
if(u .lt. 0.d0) u = u + tpi
!  Compute the radius R and velocity squared V2, and the dot
!  product RDOTV, the energy per unit mass ENERGY .
r = dsqrt(x*x + y*y + z*z)
v2 = vx*vx + vy*vy + vz*vz
v = dsqrt(v2)
vdotr = x*vx + y*vy + z*vz
energy = 0.5d0*v2 - GM/r
! ELLIPSE
  a = -0.5d0*GM/energy
  fac = 1.d0 - h2/(GM*a)
  if (fac .gt. TINY) then
     e = dsqrt ( fac ) ! eccentricity
     ec = 1.d0-r/a ! cos(eccentric anomaly)
     es = vdotr/dsqrt(GM*a) ! sin(eccentric anomaly)
     cape = datan2(es,ec) ! eccentric anomaly, E
     cw = (ec/e -e)/(1.d0 - ec)          ! cos(true anomaly)
     sw = dsqrt(1.d0-e*e)*(es/e)/(1.d0-ec) ! sin(true anomaly)
     w = datan2(sw,cw) ! true anomaly
     if(w .lt. 0.d0) w = w + tpi
     else
       e = 0.d0
       cape = u
       es = 0.d0
       w = u
     endif
     capmnq = cape - es ! mean anomaly
     peri = u - w ! arg of periapse
     peri = dmod(peri,tpi)
     if(peri .lt. 0.d0) peri = peri + tpi
return
end subroutine
