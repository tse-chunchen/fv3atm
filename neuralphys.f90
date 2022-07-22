!***************************************************************************
!
! Name: neuralphys
!
! Language: FORTRAN                           Type - MODULE
!
! Version: 1.0          Date: 03-25-20         
!     
!
! **************************************************************
!
! Module contains all subroutines require to initialize and 
! calculate ensemble of NNs for GFS model physics.
! 
! Originally created by Alex Belochitski for physics emulation
! Modified by Tse-Chun Chen for bias correction
! **************************************************************
!
       module neuralphys

        !use module_iounitdef, only : nonetf
        use sim_nc_mod,        only: open_ncfile,  &
                                     close_ncfile, &
                                     get_var2_real,&
                                     get_var1_real,&
                                     handle_err
        use module_radlw_parameters, only: sfcflw_type,  &
                                           topflw_type
        use module_radsw_parameters, only: sfcfsw_type,  &
                                           topfsw_type
        use machine,           only: kind_phys

        use omp_lib
                 
        implicit none 
#include <netcdf.inc>

        private
 
        public :: init_nn, eval_nn
                  
! Number of members in the NN ensemble

        integer, parameter     :: nn_num_of_members = 4
!
        real, parameter     :: pi = 3.1415927
! Files containing NN weights and biases

        character(*), parameter::nn_file_name(nn_num_of_members)= (/& 
             '/scratch2/BMC/gsienkf/Tse-chun.Chen/data/nn_t.nc', &
             '/scratch2/BMC/gsienkf/Tse-chun.Chen/data/nn_u.nc', &
             '/scratch2/BMC/gsienkf/Tse-chun.Chen/data/nn_v.nc', &
             '/scratch2/BMC/gsienkf/Tse-chun.Chen/data/nn_q.nc' /)
             
        !integer,      parameter::nn_ncid(nn_num_of_members)= (/101,102,103,104/)
        
        integer,    allocatable :: nn_sizes(:,:) ! nn_num_of_members x num_layers
        integer num_layers,num_layersm1

 ! Internal types and variables
        type nndata
            real(kind=4), allocatable :: w(:,:)
            real(kind=4), allocatable :: b(:)
        end type nndata

        !type nndata_1d
        !   real, allocatable :: a(:)
        !end type nndata_1d

        !type nndata_2d
        !   real, allocatable :: a(:,:)
        !end type nndata_2d
        
        !type nndata_3d
        !   real, allocatable :: a(:,:,:)
        !end type nndata_3d

! Get number of NN layers from nn_t.nc attributes
        !call open_ncfile(nn_file_name(1), 104)
        !STATUS = NF_INQ_ATTLEN (104, NF_GLOBAL, 'nn_sizes', num_layers) ! get number of NN layers
        !IF (STATUS .NE. NF_NOERR) CALL HANDLE_ERR(STATUS)
        !call close_ncfile(104)
        
        !allocate(nn_sizes(nn_num_of_members,num_layers))        

! Assume all NNs to have the same number of layers (number of neurons can be different).
        type(nndata), allocatable :: nns(:,:)

! NN hidden weights
        !type(nndata_3d) :: nn_whid(nn_num_of_members)
! NN input and output weights, hidden biases
        !type(nndata_2d) :: nn_win(nn_num_of_members),nn_wout(nn_num_of_members),nn_bhid(nn_num_of_members)
! NN output biases
        !type(nndata_1d) :: nn_bout(nn_num_of_members)

      contains
!
! Initialize TMP NNs
        
        subroutine init_nn(me) 
!
! --- This subroutine initializes NN ensemble, i.e. reads NNs coefficients
!
!   
          integer, intent(in) ::  me
          integer iin,ihid,iout,member,layer,l0,l1,STATUS,ncid !,num_layers,num_layersm1
          character(len=2) :: varname 
          
          if (me == 0) print*,'Module NEURALPHYS: Number of NNs:', &
                               nn_num_of_members

! Get number of NN layers from nn_t.nc attributes
          call open_ncfile(nn_file_name(1), ncid)
          STATUS = NF_INQ_ATTLEN (ncid, NF_GLOBAL, 'nn_sizes', num_layers) ! get number of NN layers
          IF (STATUS .NE. NF_NOERR) CALL handle_err("inq_attlen",STATUS)
          call close_ncfile(ncid)

          num_layersm1 = num_layers - 1

          allocate(nn_sizes(nn_num_of_members,num_layers))

! Assume all NNs to have the same number of layers (number of neurons can be
! different).
          allocate(nns(nn_num_of_members,num_layersm1))

! Load NNs weights and biases

          do member=1,nn_num_of_members
              call open_ncfile(nn_file_name(member), ncid)
        
              STATUS = NF_GET_ATT_INT (ncid, NF_GLOBAL, 'nn_sizes', nn_sizes(member,:))
              IF (STATUS .NE. NF_NOERR) CALL handle_err("get_att_int",STATUS)
              
              do layer=1,num_layers-1
                  l0 = nn_sizes(member,layer)
                  l1 = nn_sizes(member,layer+1)
                  allocate(nns(member,layer)%w(l0,l1))
                  allocate(nns(member,layer)%b(l1))
                  
                  write(varname, "(A1,I1)") "w",layer-1
                  call get_var2_real(ncid, varname, l0, l1, nns(member,layer)%w(:,:))
                  write(varname, "(A1,I1)") "b",layer-1
                  call get_var1_real(ncid, varname, l1, nns(member,layer)%b(:))
              end do
              call close_ncfile(ncid)
              
              if (me == 0)  print*,'Module NEURALPHYS: NN File Loaded: ', & 
                                   nn_file_name(member)
          end do
        end subroutine init_nn

! TMP emulator

        subroutine  eval_nn(                                     & 
!  Inputs:
! Prognostic variables:
     !&     pgr,phil,prsl,ugrs,vgrs,vvl,tgrs,shm,cwmr,ozmr, 	&
     &     pgr,ugrs,vgrs,tgrs,shm,  &
! Surface variables:
     !&     smc,slc,stc,tskin,canopy,hice,weasd,  &
     &     cmm,evcw,evbs,sbsno,snohf,snowc,srunoff,  &
     &     trans,tsfc,tisfc,q2m,epi,zorl,alboldxy,   &
     &     sfcflw,sfcfsw,topflw,topfsw,slmsk,        &
! Metavariables:
     !&     ftime,doy,mon,glon,glat,cosz,solcon,  &
     &     hour,doy,glon,glat,dt,  &
! Outputs:
! Prognostic variables:
     &     gu0,gv0,gt0,oshm) !ocwmr,oozmr, &
! Surface variables:
     !&     osmc,oslc,ostc,otskin,ocanopy,ohice,oweasd)


! Inputs
!          integer,  intent(in) :: jday

          real, intent(in) :: pgr,cmm,evcw,evbs,sbsno,snohf,snowc,srunoff,trans,tsfc,tisfc,q2m,epi,zorl,alboldxy,slmsk,glat,glon,hour,doy,dt
          type (sfcflw_type), intent(in) :: sfcflw
          type (sfcfsw_type), intent(in) :: sfcfsw
          type (topflw_type), intent(in) :: topflw
          type (topfsw_type), intent(in) :: topfsw

          real, intent(in):: ugrs(:),vgrs(:),tgrs(:),shm(:)
! Outputs     
          real, intent(out):: gu0(:),gv0(:),gt0(:),oshm(:)
!
! Local variables
          real  nn_input_vector(nn_sizes(1,1)),  nn_output_vector(nn_sizes(1,num_layers)) 

!             
! Create NN input vector:
!        

          nn_input_vector(1:127)  = tgrs
          nn_input_vector(128)    = log(pgr)
          nn_input_vector(129:255)= ugrs
          nn_input_vector(256:382)= vgrs
          nn_input_vector(383:509)= shm
          nn_input_vector(510)    = cmm ! or uustar; may not be the same
          nn_input_vector(511)    = evcw
          nn_input_vector(512)    = evbs
          nn_input_vector(513)    = sbsno
          nn_input_vector(514)    = snohf
          nn_input_vector(515)    = snowc !may not be the same
          nn_input_vector(516)    = srunoff !may not be the same
          nn_input_vector(517)    = trans !may not be the same
          nn_input_vector(518)    = tsfc
          nn_input_vector(519)    = tisfc
          nn_input_vector(520)    = q2m
          nn_input_vector(521)    = epi  !may not be the same
          nn_input_vector(522)    = zorl !may not be the same
          nn_input_vector(523)    = alboldxy !may not be the same
          nn_input_vector(524)    = sfcflw%dnfx0
          nn_input_vector(525)    = sfcfsw%dnfx0
          nn_input_vector(526)    = sfcflw%upfx0
          nn_input_vector(527)    = topflw%upfx0
          nn_input_vector(528)    = sfcfsw%upfx0
          nn_input_vector(529)    = topfsw%upfx0
          nn_input_vector(530)    = slmsk
          nn_input_vector(531)    = glat
          nn_input_vector(532)    = sin(glon)
          nn_input_vector(533)    = cos(glon)
          nn_input_vector(534)    = sin(2.* pi * hour/24.)
          nn_input_vector(535)    = sin(2.* pi * doy/365.)
          nn_input_vector(536)    = cos(2.* pi * hour/24.)
          nn_input_vector(537)    = cos(2.* pi * doy/365.)


          !nn_input_vector(1)  =           cos(2.* pi * doy/366.)  
          !nn_input_vector(2) =            sin(2.* pi * doy/366.)
          !nn_input_vector(3) =            cos(2.* pi * mon/12.) 
          !nn_input_vector(4) =            sin(2.* pi * mon/12.)
          !nn_input_vector(5) =            glat
          !nn_input_vector(6) =            cos(glon)
          !nn_input_vector(7) =            sin(glon)
          !nn_input_vector(8) =            cosz 
          !nn_input_vector(9:72) =        phil         ! layer geopotential height
          !nn_input_vector(73:136) =       prsl         ! layer pressure
          !nn_input_vector(137) =          pgr          ! surface pressure 
          !nn_input_vector(138:201) =      ugrs         ! u component of layer wind
          !nn_input_vector(202:265) =      vgrs         ! v component of layer wind
          !nn_input_vector(266:329) =      vvl          ! layer mean vertical velocity
          !nn_input_vector(330:393) =      tgrs         ! layer mean temperature
          !nn_input_vector(394:430) =      shm(1:37)    ! specific humidity
          !nn_input_vector(431:473) =      cwmr(1:43)   ! cloud water mixing ratio
          !nn_input_vector(474:505) =      ozmr(33:64)  ! ozone mixing ratio
          !nn_input_vector(506) =          solcon       ! solar constant

!             
! Call NN computation                     
          call compute_nn(nn_input_vector,nn_output_vector) !,nn_num_of_members,& 
                  !nn_w1,nn_w2,nn_b1,nn_b2,nn_hid)
!
! Unpack NN output vector
          gu0 	       = nn_output_vector(128:254)*dt + ugrs   ! u component of layer wind
          gv0          = nn_output_vector(255:381)*dt + vgrs ! v component of layer wind
          gt0          = nn_output_vector(1:127)*dt + tgrs! layer mean temperature 
!          where (nn_output_vector(193:304) < 0.) nn_output_vector(193:304) = 0.
          !oshm(1:37)   = nn_output_vector(193:229)
          !oshm(38:64)  = 0.
          oshm         = nn_output_vector(382:508)*dt + shm                      ! specific humidity
          !ocwmr(1:43)  = nn_output_vector(230:272) 
          !ocwmr(44:64) = 0. 
          !ocwmr        = ocwmr + cwmr                    ! cloud water mixing ratio
          !oozmr(33:64) = nn_output_vector(273:304)    
          !oozmr(1:32)  = 0.
          !oozmr        = oozmr + ozmr                    ! ozone mixing ratio
!
        end subroutine eval_nn
        
        subroutine  compute_nn(X,Y) !,num_of_members,w1,w2,b1,b2,nhid)

      
 !  Input:
 !            X(IN) NN input vector 
          !integer, intent(in) :: num_of_members,nhid(num_of_members)
          real, intent(in)::X(:)
          !type(nndata_2d), intent(in) :: w1(num_of_members),w2(num_of_members)
          !type(nndata_1d), intent(in) :: b1(num_of_members),b2(num_of_members)
         

 !   Ouput:
 !            Y(OUT) NN output vector (composition coefficients for SNN)

          real, intent(out):: Y(:)

! Local variables 
          integer i, nout
          real, allocatable :: x_tmp1(:), x_tmp2(:)! x2(:),x3(:)
          integer member,layer, l0, l1
          
          do member = 1,nn_num_of_members              
              do layer = 1,num_layers-2 ! loop from 1st to final-1 layer
                  l0 = nn_sizes(member,layer)
                  l1 = nn_sizes(member,layer+1)
                  allocate(x_tmp2(l1))
                  
                  if (layer == 1) then
                      allocate(x_tmp1(l0))
                      x_tmp1=X
                  endif
! Internal layers                  
!$OMP PARALLEL default (shared) private (i)                                                                           
!$OMP DO
                  do i = 1,l1
                      x_tmp2(i) = tanh(sum(x_tmp1*nns(member,layer)%w(:,i))+nns(member,layer)%b(i))
                  end do
!$OMP END DO                                                                                                          
!$OMP END PARALLEL 
                  
                  deallocate(x_tmp1)
                  allocate(x_tmp1(l1))
                  x_tmp1 = x_tmp2
                  deallocate(x_tmp2)
              end do
              
              !l0 = nn_sizes(member,num_layers-1)
              !l1 = nn_sizes(member,num_layers)
              !allocate(x_tmp2(l1))
              
! Output layer              
!$OMP PARALLEL default (shared) private (i)                                                                           
!$OMP DO              
              do i = 1,127 !l1
                  Y(member*127+i) = sum(x_tmp1*nns(member,layer)%w(:,i))+nns(member,layer)%b(i)
              end do
!$OMP END DO                                                                                                          
!$OMP END PARALLEL 
              
              deallocate(x_tmp1)
          end do    
                  

          !nout=size(Y)

          !Y = 0.

          !allocate(x3(nout)) 

          !do member = 1,num_of_members
  
             !allocate(x2(nhid(member))) 

! Calculate neurons in the hidden layer

!             forall(i = 1:nhid(member)) x2(i)= tanh(sum(X*w1(member)%a(:,i))+  & 
!                                               b1(member)%a(i))

!!$OMP PARALLEL default (shared) private (i)                                                                           
!!$OMP DO                                                                                                              
             !do i = 1,nhid(member)

             !   x2(i)= tanh(sum(X*w1(member)%a(:,i)) + b1(member)%a(i))

             !enddo
!!$OMP END DO                                                                                                          
!!$OMP END PARALLEL 

! Calculate NN output 

!             forall(i=1:nout) x3(i)= sum(w2(member)%a(:,i)*x2) + b2(member)%a(i)

!!$OMP PARALLEL default (shared) private (i)                                                                           
!!$OMP DO                                                                                                              
             !do i=1,nout

             !   x3(i)= sum(w2(member)%a(:,i)*x2) + b2(member)%a(i)

             !enddo
!!$OMP END DO                                                                                                          
!!$OMP END PARALLEL 

             !Y = Y + x3 
             
             !deallocate(x2)
         
          !end do                    ! member

          !deallocate(x3)
      
          !Y = Y / num_of_members
         
    end  subroutine  compute_nn

  end module neuralphys
