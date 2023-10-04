      subroutine TimeMarcher
      use param
      use local_arrays
      use mpih
      use mpi_param, only: kstart,kend
      use mls_param
      use local_aux
      use mls_local, only: for_xc, for_yc, for_zc
      use stat_arrays, only: vxvyvz_rms_vol
      use pointparticle
      implicit none
      integer :: ns, inp, ntr, nsub,i
      real :: radius,zpost,sumf
      

!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     MLS Routines
! 

      nsstep = 1

!     Computing normals and triangle properties
      do inp=1,Nparticle

       do i = 1,maxnv
         radius=sqrt((xyz(2,i,inp)-0.5)**2.+(xyz(3,i,inp)-0.5)**2.)
        
         xyzv(1,i,inp)=wirev
         xyzv(2,i,inp)=0.
         xyzv(3,i,inp)=0.
        
         xyz(1:3,i,inp)=xyz(1:3,i,inp)+(dt/float(nsstep)) &
                             *(xyzv(1:3,i,inp))

       end do

        call convert_geo(maxnv,maxne,maxnf,xyz(:,:,inp),xyzv(:,:,inp), &
                  xyza(:,:,inp),vert_of_face(:,:,inp), &
                  tri_ver(:,:,inp),tri_bar(:,:,inp), &
                  vel_tri(:,:,inp),acc_tri(:,:,inp))

        call calculate_area (Surface(inp),maxnv,maxnf,xyz(:,:,inp), &
                              vert_of_face(:,:,inp),sur(:,inp))
      end do

!     Find indices of all centroids of particles and bounding box
      call findindices

!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Point particle routines

      if(withppart) then 
       if(time.ge.timeONp.and.ONparticle.eq.0)then
         call InitPointP
         for_x_part=0.d0; for_y_part=0.d0; for_z_part=0.d0
       end if
       call UpdatePointParticlePosition
       call storeold
       if(time.ge.timeONp) ONparticle=1
      end if
!    
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!
!   TIME INTEGRATION : implicit viscous
!                                                                       

      beta=dt/ren*0.5d0

      do ns=1,nsst                                                 
        al=alm(ns)
        ga=gam(ns)
        ro=rom(ns)

!     ============================================
!     Start calculation of intermediate, non-solenoidal velocity
!     ============================================

        call ExplicitTermsVX
        call ExplicitTermsVY
        call ExplicitTermsVZ

        call ImplicitAndUpdateVX
        call ImplicitAndUpdateVY     
        call ImplicitAndUpdateVZ     

        call update_both_ghosts(n1,n2,vx,kstart,kend)
        call update_both_ghosts(n1,n2,vy,kstart,kend)
        call update_both_ghosts(n1,n2,vz,kstart,kend)
 

!     ============================================
!     Forcing for fluid-structure interaction
!     ============================================

       if(mlsforcing) then
        for_xc=0.d0 ; for_yc=0.d0 ; for_zc=0.d0
        call CalcMLSForce
        call MLSForceVelocity
       end if

!     ======================================================
!     End of MLS forcing & intermediate velocity field calculation
!     ======================================================

        call update_upper_ghost(n1,n2,vz)

!     ======================================================
!     Start pressure correction step
!     ======================================================

        call CalcLocalDivergence !here only vz ghost(up) cell are needed
        call SolvePressureCorrection

        call update_both_ghosts(n1,n2+1,dph,kstart,kend)
        
        call CorrectVelocity                 !! SOLENOIDAL VEL FIELD
        call CorrectPressure                         !! PRESSURE FIELD

        call update_both_ghosts(n1,n2,vx,kstart,kend)
        call update_both_ghosts(n1,n2,vy,kstart,kend)
        call update_both_ghosts(n1,n2,vz,kstart,kend)
        call update_both_ghosts(n1,n2,pr,kstart,kend)

!     ======================================================
!     End pressure correction
!     ======================================================

        enddo


      if(withppart) then

!CS   Calculate material derivative at staggered location
      call CalcMaterialDerivative
      call update_both_ghosts(n1,n2,matderx,kstart,kend)
      call update_both_ghosts(n1,n2,matdery,kstart,kend)
      call update_both_ghosts(n1,n2,matderz,kstart,kend)

!CS   Center velocity and material derivative
      call CentreVariables

      call CalcVorticity
      call update_upper_ghost(n1m,n2m,vorx)
      call update_upper_ghost(n1m,n2m,vory)
      call update_upper_ghost(n1m,n2m,vorz)

!CS   Now solve for velocity of particles
      if(ONparticle.eq.1) call CalcPointPVel

!CS   Now reset forcing
       for_x_part(:,:,:)=0.d0
       for_y_part(:,:,:)=0.d0
       for_z_part(:,:,:)=0.d0
      end if



!m================================       
!m    CALLS FOR STRUCTURAL SOLVER

    
      if(solvestructure)then

       do nsub = 1,nsstep

          do inp=1,Nparticle

!     Computing normals and triangle properties
        call calculate_normal(maxnv,maxnf, &
                  xyz(:,:,inp),vert_of_face(:,:,inp), &
                  tri_nor(:,:,inp))

          end do
      
!     Compute pressure and viscous loads
       call mlsStruc_closed_rigid

!     Reduce the forces from all processors over each particle
       do ntr=1,maxnf
        call mpi_globalsum_double_arr(press_face(ntr,:),Nparticle)
        call mpi_globalsum_double_arr(vforc_face(ntr,:),Nparticle)
        call mpi_globalsum_double_arr(sca_nod1(ntr,:),Nparticle)
        call mpi_globalsum_double_arr(sca_nod2(ntr,:),Nparticle)
        call mpi_globalsum_double_arr(sca_nod3(ntr,:),Nparticle)
       end do

       call MpiBarrier

!     ------ Output simple statistics if needed ------------
       if(ismaster)then
        if(maxval(xyzv(:,:,:)).gt.10.0)then
         print*,'Particle exploded', maxval(xyzv(:,:,:))
         stop
         call MpiAbort
        end if
       end if


       if(ismaster) then
        sumf=0.
        do inp=1,Nparticle
        do ntr=1,maxnf
         zpost=tri_bar(3,ntr,inp)
         if((zpost.gt.0.3).and.(zpost.lt.2.7)) then
          sumf=sumf+press_face(ntr,1)
         end if
        end do
        end do

        write(113,*)time,sumf

!        write(114,*)time,sum(press_face(:,1))/vxvyvz_rms_vol**2

!  115 is viscf.out, added for use with closed_rigid objects
!        write(115,*)time,sum(vforc_face(:,1))
       end if
!     --------------------------------------------------------------

       end do
      end if


!    END OF  STRUCTURAL SOLVER
!================================

 
 

      return                                                            
      end                                                               
