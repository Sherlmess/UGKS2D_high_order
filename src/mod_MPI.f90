module mpi_wrapper
    use mpi
    use mpi_para
    use flow_data
    implicit none
    
contains
      
    !partion physical space
    subroutine partion_mpi()
        implicit none
        integer(KINT)           :: ip,ipx1,ipx2,ipy1,ipy2,ki,kj,np_size

        !total processer
        NP_XY = NP_X*NP_Y
        call MPI_Comm_Size(MPI_COMM_WORLD,np_size,ierr)      ! get total np size
        if(np_size .ne. NP_XY) then
            if(myID .eq. 0) then
                print*, 'Number of Total Procs is ',np_size
                print*, 'Number of Total procs is not equal to NP_X*NP_Y!!!'
            endif
            stop
        endif
        
        !current processer
        IP_XY = myID
        IP_X  = mod(IP_XY,NP_X)                   ! block index in x-direction 
        IP_Y  = IP_XY/NP_X                        ! block index in y-direction
        
        !allocate array
        allocate(Block_XY(0:NP_XY-1))
        
        !block size and index
        Nx = Nx_all/NP_X
        Ny = Ny_all/NP_Y
        Block_XY(:)%NI = Nx
        Block_XY(:)%NJ = Ny
        do kj=0,NP_Y-1                                  !consideration for not-divisible
        !load balance: the first mod(Nx_all,NP_X) blocks in x-direction have one more cell, and the first mod(Ny_all,NP_Y) blocks in y-direction have one more cell
            do ki=0,NP_X-1                              
                ip = kj*NP_X+ki
                if(ki<mod(Nx_all,NP_X)) Block_XY(ip)%NI=Block_XY(ip)%NI+1
                if(kj<mod(Ny_all,NP_Y)) Block_XY(ip)%NJ=Block_XY(ip)%NJ+1
                Block_XY(ip)%Origin_I = int(Nx_all/NP_X)*ki+min(ki,mod(Nx_all,NP_X))+1
                Block_XY(ip)%Origin_J = int(Ny_all/NP_Y)*kj+min(kj,mod(Ny_all,NP_Y))+1
            enddo
        enddo
        
        !cell numbers of current block
        Nx = Block_XY(IP_XY)%NI
        Ny = Block_XY(IP_XY)%NJ
        
        !communication node with current block
        ipx1        = mymod(IP_X-1,NP_X)                !ID_X_Minus:  -1 proc in x-direction
        ipx2        = mymod(IP_X+1,NP_X)                !ID_X_Plus:   +1 proc in x-direction
        ID_X_Minus  = IP_Y*NP_X+ipx1                    
        ID_X_Plus   = IP_Y*NP_X+ipx2                    
        if(bc_type(1) /= BCconst_periodic .and. IP_X == 0) then         !If IP_X=0: NP_X-1 node send mesg to 0 node for periodic boundary, else no comm.
            ID_X_Minus = MPI_PROC_NULL
        endif
        if(bc_type(2) /= BCconst_periodic .and. IP_X == NP_X-1) then    !If IP_X=NP_X-1: 0 node send mesg to NP_X-1 node for periodic boundary, else no comm.
            ID_X_Plus  = MPI_PROC_NULL  
        endif
      
        ipy1        = mymod(IP_Y-1,NP_Y)                !ID_Y_Minus:  -1 proc in y-direction
        ipy2        = mymod(IP_Y+1,NP_Y)                !ID_Y_Plus:   +1 proc in y-direction
        ID_Y_Minus  = ipy1*NP_X+IP_X                    
        ID_Y_Plus   = ipy2*NP_X+IP_X                    
        if(bc_type(3) /= BCconst_periodic .and. IP_Y == 0) then         !If IP_Y=0: NP_Y-1 node send mesg to 0 node for periodic boundary, else no comm.
            ID_Y_Minus = MPI_PROC_NULL
        endif
        if(bc_type(4) /= BCconst_periodic .and. IP_Y == NP_Y-1) then    !If IP_Y=NP_Y-1: 0 node send mesg to NP_Y-1 node for periodic boundary, else no comm.
            ID_Y_Plus  = MPI_PROC_NULL
        endif
        
        contains
    
        function mymod(i,n)
            implicit none
            integer mymod,i,n
            if(i.lt.0) then
                mymod=i+n
            else if (i.gt.n-1) then
                mymod=i-n
            else
                mymod=i
            endif
        end function mymod
        
    end subroutine partion_mpi


    subroutine MPIUpdate_CellCoordinates()
    
        implicit none
        integer                     :: i,j,k1,Ndata_send
        real(KREAL),allocatable     :: tmp_send1(:,:),tmp_send2(:,:),tmp_recv1(:,:),tmp_recv2(:,:)

        !Sendrecv: ctr%x, y, length(1), length(2), area 
        Ndata_send=5
        
        !------------------------------------------
        !>I direction
        !------------------------------------------
        allocate(tmp_send1(LAP*ny,Ndata_send),tmp_send2(LAP*ny,Ndata_send),tmp_recv1(LAP*ny,Ndata_send),tmp_recv2(LAP*ny,Ndata_send))
        do j=1,ny
            do i=1,LAP
                k1=(j-1)*LAP+i
                tmp_send1(k1,1)=ctr(i,j)%x
                tmp_send1(k1,2)=ctr(i,j)%y
                tmp_send1(k1,3)=ctr(i,j)%length(1)
                tmp_send1(k1,4)=ctr(i,j)%length(2)
                tmp_send1(k1,5)=ctr(i,j)%area
                
                tmp_send2(k1,1)=ctr(nx-LAP+i,j)%x
                tmp_send2(k1,2)=ctr(nx-LAP+i,j)%y
                tmp_send2(k1,3)=ctr(nx-LAP+i,j)%length(1)
                tmp_send2(k1,4)=ctr(nx-LAP+i,j)%length(2)
                tmp_send2(k1,5)=ctr(nx-LAP+i,j)%area
            enddo
        enddo

        call MPI_Sendrecv(tmp_send1, LAP*ny*Ndata_send, MPI_KREAL, ID_X_minus, 7000, &
                          tmp_recv2, LAP*ny*Ndata_send, MPI_KREAL, ID_X_plus,  7000, MPI_COMM_WORLD, Status, ierr)
        call MPI_Sendrecv(tmp_send2, LAP*ny*Ndata_send, MPI_KREAL, ID_X_plus,  6000,   &
                          tmp_recv1, LAP*ny*Ndata_send, MPI_KREAL, ID_X_minus, 6000, MPI_COMM_WORLD, Status, ierr)

        if(ID_X_minus .ne. MPI_PROC_NULL) then
            do j=1,ny
                do i=1,LAP
                    k1=(j-1)*LAP+i
                    ctr(i-LAP,j)%x          = tmp_recv1(k1,1)
                    ctr(i-LAP,j)%y          = tmp_recv1(k1,2)
                    ctr(i-LAP,j)%length(1)  = tmp_recv1(k1,3)
                    ctr(i-LAP,j)%length(2)  = tmp_recv1(k1,4)
                    ctr(i-LAP,j)%area       = tmp_recv1(k1,5)
                enddo
            enddo
        endif

        if(ID_X_plus .ne. MPI_PROC_NULL) then
            do j=1,ny
                do i=1,LAP
                    k1=(j-1)*LAP+i
                    ctr(nx+i,j)%x           = tmp_recv2(k1,1)
                    ctr(nx+i,j)%y           = tmp_recv2(k1,2)
                    ctr(nx+i,j)%length(1)   = tmp_recv2(k1,3)
                    ctr(nx+i,j)%length(2)   = tmp_recv2(k1,4)
                    ctr(nx+i,j)%area        = tmp_recv2(k1,5)
                enddo
            enddo
        endif
        deallocate(tmp_send1,tmp_send2,tmp_recv1,tmp_recv2)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
        
        !------------------------------------------
        !>J direction
        !------------------------------------------
        allocate(tmp_send1(LAP*nx,Ndata_send),tmp_send2(LAP*nx,Ndata_send),tmp_recv1(LAP*nx,Ndata_send),tmp_recv2(LAP*nx,Ndata_send))
        do j=1,LAP
            do i=1,nx
                k1=(j-1)*nx+i
                tmp_send1(k1,1)=ctr(i,j)%x
                tmp_send1(k1,2)=ctr(i,j)%y
                tmp_send1(k1,3)=ctr(i,j)%length(1)
                tmp_send1(k1,4)=ctr(i,j)%length(2)
                tmp_send1(k1,5)=ctr(i,j)%area
                
                tmp_send2(k1,1)=ctr(i,ny+j-LAP)%x
                tmp_send2(k1,2)=ctr(i,ny+j-LAP)%y
                tmp_send2(k1,3)=ctr(i,ny+j-LAP)%length(1)
                tmp_send2(k1,4)=ctr(i,ny+j-LAP)%length(2)
                tmp_send2(k1,5)=ctr(i,ny+j-LAP)%area
            enddo
        enddo

        call MPI_Sendrecv(tmp_send1, LAP*nx*Ndata_send, MPI_KREAL, ID_Y_minus, 9000, &
                          tmp_recv2, LAP*nx*Ndata_send, MPI_KREAL, ID_Y_plus,  9000, MPI_COMM_WORLD, Status, ierr)
        call MPI_Sendrecv(tmp_send2, LAP*nx*Ndata_send, MPI_KREAL, ID_Y_plus,  8000,   &
                          tmp_recv1, LAP*nx*Ndata_send, MPI_KREAL, ID_Y_minus, 8000, MPI_COMM_WORLD, Status, ierr)


        if(ID_Y_minus .ne. MPI_PROC_NULL) then
            do j=1,LAP
                do i=1,nx
                    k1=(j-1)*nx+i
                    ctr(i,j-LAP)%x          = tmp_recv1(k1,1)
                    ctr(i,j-LAP)%y          = tmp_recv1(k1,2)
                    ctr(i,j-LAP)%length(1)  = tmp_recv1(k1,3)
                    ctr(i,j-LAP)%length(2)  = tmp_recv1(k1,4)
                    ctr(i,j-LAP)%area       = tmp_recv1(k1,5)
                enddo
            enddo
        endif

        if(ID_Y_plus .ne. MPI_PROC_NULL) then
            do j=1,LAP
                do i=1,nx
                    k1=(j-1)*nx+i
                    ctr(i,ny+j)%x           = tmp_recv2(k1,1)
                    ctr(i,ny+j)%y           = tmp_recv2(k1,2)
                    ctr(i,ny+j)%length(1)   = tmp_recv2(k1,3)
                    ctr(i,ny+j)%length(2)   = tmp_recv2(k1,4)
                    ctr(i,ny+j)%area        = tmp_recv2(k1,5)
                enddo
            enddo
        endif
        deallocate(tmp_send1,tmp_send2,tmp_recv1,tmp_recv2)
        
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    end subroutine MPIUpdate_CellCoordinates
    
    
    subroutine MPIUpdate_CellConVars()
    
        implicit none
        integer                     :: i,j,k1,Ndata_send
        real(KREAL),allocatable     :: tmp_send1(:,:),tmp_send2(:,:),tmp_recv1(:,:),tmp_recv2(:,:)

        !Sendrecv: ctr%ConVars(1:4) 
        Ndata_send=4
        
        !------------------------------------------
        !>I direction
        !------------------------------------------
        allocate(tmp_send1(LAP*ny,Ndata_send),tmp_send2(LAP*ny,Ndata_send),tmp_recv1(LAP*ny,Ndata_send),tmp_recv2(LAP*ny,Ndata_send))
        do j=1,ny
            do i=1,LAP
                k1=(j-1)*LAP+i
                tmp_send1(k1,1:4)=ctr(i,j)%ConVars(1:4)
                tmp_send2(k1,1:4)=ctr(nx-LAP+i,j)%ConVars(1:4)
            enddo
        enddo

        call MPI_Sendrecv(tmp_send1, LAP*ny*Ndata_send, MPI_KREAL, ID_X_minus, 7000, &
                          tmp_recv2, LAP*ny*Ndata_send, MPI_KREAL, ID_X_plus,  7000, MPI_COMM_WORLD, Status, ierr)
        call MPI_Sendrecv(tmp_send2, LAP*ny*Ndata_send, MPI_KREAL, ID_X_plus,  6000,   &
                          tmp_recv1, LAP*ny*Ndata_send, MPI_KREAL, ID_X_minus, 6000, MPI_COMM_WORLD, Status, ierr)

        if(ID_X_minus .ne. MPI_PROC_NULL) then
            do j=1,ny
                do i=1,LAP
                    k1=(j-1)*LAP+i
                    ctr(i-LAP,j)%ConVars(1:4) = tmp_recv1(k1,1:4)
                enddo
            enddo
        endif

        if(ID_X_plus .ne. MPI_PROC_NULL) then
            do j=1,ny
                do i=1,LAP
                    k1=(j-1)*LAP+i
                    ctr(nx+i,j)%ConVars(1:4) = tmp_recv2(k1,1:4)
                enddo
            enddo
        endif
        deallocate(tmp_send1,tmp_send2,tmp_recv1,tmp_recv2)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
        
        !------------------------------------------
        !>J direction
        !------------------------------------------
        allocate(tmp_send1(LAP*nx,Ndata_send),tmp_send2(LAP*nx,Ndata_send),tmp_recv1(LAP*nx,Ndata_send),tmp_recv2(LAP*nx,Ndata_send))
        do j=1,LAP
            do i=1,nx
                k1=(j-1)*nx+i
                tmp_send1(k1,1:4)=ctr(i,j)%ConVars(1:4)
                tmp_send2(k1,1:4)=ctr(i,ny+j-LAP)%ConVars(1:4)
            enddo
        enddo

        call MPI_Sendrecv(tmp_send1, LAP*nx*Ndata_send, MPI_KREAL, ID_Y_minus, 9000, &
                          tmp_recv2, LAP*nx*Ndata_send, MPI_KREAL, ID_Y_plus,  9000, MPI_COMM_WORLD, Status, ierr)
        call MPI_Sendrecv(tmp_send2, LAP*nx*Ndata_send, MPI_KREAL, ID_Y_plus,  8000,   &
                          tmp_recv1, LAP*nx*Ndata_send, MPI_KREAL, ID_Y_minus, 8000, MPI_COMM_WORLD, Status, ierr)


        if(ID_Y_minus .ne. MPI_PROC_NULL) then
            do j=1,LAP
                do i=1,nx
                    k1=(j-1)*nx+i
                    ctr(i,j-LAP)%ConVars(1:4) = tmp_recv1(k1,1:4)
                enddo
            enddo
        endif

        if(ID_Y_plus .ne. MPI_PROC_NULL) then
            do j=1,LAP
                do i=1,nx
                    k1=(j-1)*nx+i
                    ctr(i,ny+j)%ConVars(1:4) = tmp_recv2(k1,1:4)
                enddo
            enddo
        endif
        deallocate(tmp_send1,tmp_send2,tmp_recv1,tmp_recv2)
        
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    end subroutine MPIUpdate_CellConVars
    
    
    subroutine MPIUpdate_CellVDF()
        implicit none
        integer                                     :: i,j,k1,Ndata_send
        real(KREAL),allocatable,dimension(:,:,:)    :: tmp_send1,tmp_send2,tmp_recv1,tmp_recv2

        !Sendrecv: flattened cell-center VDFs h_vec/b_vec
        Ndata_send=2
        
        !------------------------------------------
        !>I direction
        !------------------------------------------
        allocate(tmp_send1(LAP*ny,Nvel,Ndata_send),tmp_send2(LAP*ny,Nvel,Ndata_send),&
                 tmp_recv1(LAP*ny,Nvel,Ndata_send),tmp_recv2(LAP*ny,Nvel,Ndata_send))
        do j=1,ny
            do i=1,LAP
                k1=(j-1)*LAP+i
                tmp_send1(k1,1:Nvel,1)=ctr(i,j)%h_vec(1:Nvel)
                tmp_send1(k1,1:Nvel,2)=ctr(i,j)%b_vec(1:Nvel)
                tmp_send2(k1,1:Nvel,1)=ctr(nx-LAP+i,j)%h_vec(1:Nvel)
                tmp_send2(k1,1:Nvel,2)=ctr(nx-LAP+i,j)%b_vec(1:Nvel)
            enddo
        enddo

        call MPI_Sendrecv(tmp_send1, LAP*ny*Nvel*Ndata_send, MPI_KREAL, ID_X_minus, 7000, &
                          tmp_recv2, LAP*ny*Nvel*Ndata_send, MPI_KREAL, ID_X_plus,  7000, MPI_COMM_WORLD, Status, ierr)
        call MPI_Sendrecv(tmp_send2, LAP*ny*Nvel*Ndata_send, MPI_KREAL, ID_X_plus,  6000,   &
                          tmp_recv1, LAP*ny*Nvel*Ndata_send, MPI_KREAL, ID_X_minus, 6000, MPI_COMM_WORLD, Status, ierr)

        if(ID_X_minus .ne. MPI_PROC_NULL) then
            do j=1,ny
                do i=1,LAP
                    k1=(j-1)*LAP+i
                    ctr(i-LAP,j)%h_vec(1:Nvel) = tmp_recv1(k1,1:Nvel,1)
                    ctr(i-LAP,j)%b_vec(1:Nvel) = tmp_recv1(k1,1:Nvel,2)
                enddo
            enddo
        endif

        if(ID_X_plus .ne. MPI_PROC_NULL) then
            do j=1,ny
                do i=1,LAP
                    k1=(j-1)*LAP+i
                    ctr(nx+i,j)%h_vec(1:Nvel) = tmp_recv2(k1,1:Nvel,1)
                    ctr(nx+i,j)%b_vec(1:Nvel) = tmp_recv2(k1,1:Nvel,2)
                enddo
            enddo
        endif
        deallocate(tmp_send1,tmp_send2,tmp_recv1,tmp_recv2)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
        
        !------------------------------------------
        !>J direction
        !------------------------------------------
        allocate(tmp_send1(LAP*nx,Nvel,Ndata_send),tmp_send2(LAP*nx,Nvel,Ndata_send),&
                 tmp_recv1(LAP*nx,Nvel,Ndata_send),tmp_recv2(LAP*nx,Nvel,Ndata_send))
        do j=1,LAP
            do i=1,nx
                k1=(j-1)*nx+i
                tmp_send1(k1,1:Nvel,1)=ctr(i,j)%h_vec(1:Nvel)
                tmp_send1(k1,1:Nvel,2)=ctr(i,j)%b_vec(1:Nvel)
                tmp_send2(k1,1:Nvel,1)=ctr(i,ny+j-LAP)%h_vec(1:Nvel)
                tmp_send2(k1,1:Nvel,2)=ctr(i,ny+j-LAP)%b_vec(1:Nvel)
            enddo
        enddo

        call MPI_Sendrecv(tmp_send1, LAP*nx*Nvel*Ndata_send, MPI_KREAL, ID_Y_minus, 9000, &
                          tmp_recv2, LAP*nx*Nvel*Ndata_send, MPI_KREAL, ID_Y_plus,  9000, MPI_COMM_WORLD, Status, ierr)
        call MPI_Sendrecv(tmp_send2, LAP*nx*Nvel*Ndata_send, MPI_KREAL, ID_Y_plus,  8000,   &
                          tmp_recv1, LAP*nx*Nvel*Ndata_send, MPI_KREAL, ID_Y_minus, 8000, MPI_COMM_WORLD, Status, ierr)


        if(ID_Y_minus .ne. MPI_PROC_NULL) then
            do j=1,LAP
                do i=1,nx
                    k1=(j-1)*nx+i
                    ctr(i,j-LAP)%h_vec(1:Nvel) = tmp_recv1(k1,1:Nvel,1)
                    ctr(i,j-LAP)%b_vec(1:Nvel) = tmp_recv1(k1,1:Nvel,2)
                enddo
            enddo
        endif

        if(ID_Y_plus .ne. MPI_PROC_NULL) then
            do j=1,LAP
                do i=1,nx
                    k1=(j-1)*nx+i
                    ctr(i,ny+j)%h_vec(1:Nvel) = tmp_recv2(k1,1:Nvel,1)
                    ctr(i,ny+j)%b_vec(1:Nvel) = tmp_recv2(k1,1:Nvel,2)
                enddo
            enddo
        endif
        deallocate(tmp_send1,tmp_send2,tmp_recv1,tmp_recv2)
        
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    end subroutine MPIUpdate_CellVDF
    
    
    !----------------------------------------------------------------------------
    !>Updata the information of conner points: for MultiDimension Flux computation
    !-----------------------------------------------------------------------------
    subroutine MPIUpdate_Conner_Coordinates()
    
        implicit none
        integer                     :: i,j,Ndata_send
        real(KREAL),allocatable     :: tmp_send1(:,:,:),tmp_send2(:,:,:),tmp_recv1(:,:,:),tmp_recv2(:,:,:)

        !Sendrecv: ctr%x, y, length(1), length(2), area 
        Ndata_send=5
        
        !------------------------------------------
        !>I direction
        !------------------------------------------
        allocate(tmp_send1(LAP,LAP*2,Ndata_send),tmp_send2(LAP,LAP*2,Ndata_send))
        allocate(tmp_recv1(LAP,LAP*2,Ndata_send),tmp_recv2(LAP,LAP*2,Ndata_send))
        
        do i=1,LAP
            tmp_send1(i,1:LAP,1)=ctr(i,1-LAP:0)%x
            tmp_send1(i,1:LAP,2)=ctr(i,1-LAP:0)%y
            tmp_send1(i,1:LAP,3)=ctr(i,1-LAP:0)%length(1)
            tmp_send1(i,1:LAP,4)=ctr(i,1-LAP:0)%length(2)
            tmp_send1(i,1:LAP,5)=ctr(i,1-LAP:0)%area
            tmp_send1(i,LAP+1:LAP+LAP,1)=ctr(i,ny+1:ny+LAP)%x
            tmp_send1(i,LAP+1:LAP+LAP,2)=ctr(i,ny+1:ny+LAP)%y
            tmp_send1(i,LAP+1:LAP+LAP,3)=ctr(i,ny+1:ny+LAP)%length(1)
            tmp_send1(i,LAP+1:LAP+LAP,4)=ctr(i,ny+1:ny+LAP)%length(2)
            tmp_send1(i,LAP+1:LAP+LAP,5)=ctr(i,ny+1:ny+LAP)%area
                
            tmp_send2(i,1:LAP,1)=ctr(nx-LAP+i,1-LAP:0)%x
            tmp_send2(i,1:LAP,2)=ctr(nx-LAP+i,1-LAP:0)%y
            tmp_send2(i,1:LAP,3)=ctr(nx-LAP+i,1-LAP:0)%length(1)
            tmp_send2(i,1:LAP,4)=ctr(nx-LAP+i,1-LAP:0)%length(2)
            tmp_send2(i,1:LAP,5)=ctr(nx-LAP+i,1-LAP:0)%area
            tmp_send2(i,LAP+1:LAP+LAP,1)=ctr(nx-LAP+i,ny+1:ny+LAP)%x
            tmp_send2(i,LAP+1:LAP+LAP,2)=ctr(nx-LAP+i,ny+1:ny+LAP)%y
            tmp_send2(i,LAP+1:LAP+LAP,3)=ctr(nx-LAP+i,ny+1:ny+LAP)%length(1)
            tmp_send2(i,LAP+1:LAP+LAP,4)=ctr(nx-LAP+i,ny+1:ny+LAP)%length(2)
            tmp_send2(i,LAP+1:LAP+LAP,5)=ctr(nx-LAP+i,ny+1:ny+LAP)%area
        enddo

        call MPI_Sendrecv(tmp_send1, LAP*LAP*2*Ndata_send, MPI_KREAL, ID_X_minus, 7000, &
                          tmp_recv2, LAP*LAP*2*Ndata_send, MPI_KREAL, ID_X_plus,  7000, MPI_COMM_WORLD, Status, ierr)
        call MPI_Sendrecv(tmp_send2, LAP*LAP*2*Ndata_send, MPI_KREAL, ID_X_plus,  6000,   &
                          tmp_recv1, LAP*LAP*2*Ndata_send, MPI_KREAL, ID_X_minus, 6000, MPI_COMM_WORLD, Status, ierr)

        if(ID_X_minus .ne. MPI_PROC_NULL) then
            do i=1,LAP
                ctr(i-LAP,1-LAP:0)%x          = tmp_recv1(i,1:LAP,1)
                ctr(i-LAP,1-LAP:0)%y          = tmp_recv1(i,1:LAP,2)
                ctr(i-LAP,1-LAP:0)%length(1)  = tmp_recv1(i,1:LAP,3)
                ctr(i-LAP,1-LAP:0)%length(2)  = tmp_recv1(i,1:LAP,4)
                ctr(i-LAP,1-LAP:0)%area       = tmp_recv1(i,1:LAP,5)
                ctr(i-LAP,ny+1:ny+LAP)%x          = tmp_recv1(i,LAP+1:LAP+LAP,1)
                ctr(i-LAP,ny+1:ny+LAP)%y          = tmp_recv1(i,LAP+1:LAP+LAP,2)
                ctr(i-LAP,ny+1:ny+LAP)%length(1)  = tmp_recv1(i,LAP+1:LAP+LAP,3)
                ctr(i-LAP,ny+1:ny+LAP)%length(2)  = tmp_recv1(i,LAP+1:LAP+LAP,4)
                ctr(i-LAP,ny+1:ny+LAP)%area       = tmp_recv1(i,LAP+1:LAP+LAP,5)
            enddo
        endif

        if(ID_X_plus .ne. MPI_PROC_NULL) then
            do i=1,LAP
                ctr(nx+i,1-LAP:0)%x           = tmp_recv2(i,1:LAP,1)
                ctr(nx+i,1-LAP:0)%y           = tmp_recv2(i,1:LAP,2)
                ctr(nx+i,1-LAP:0)%length(1)   = tmp_recv2(i,1:LAP,3)
                ctr(nx+i,1-LAP:0)%length(2)   = tmp_recv2(i,1:LAP,4)
                ctr(nx+i,1-LAP:0)%area        = tmp_recv2(i,1:LAP,5)
                ctr(nx+i,ny+1:ny+LAP)%x           = tmp_recv2(i,LAP+1:LAP+LAP,1)
                ctr(nx+i,ny+1:ny+LAP)%y           = tmp_recv2(i,LAP+1:LAP+LAP,2)
                ctr(nx+i,ny+1:ny+LAP)%length(1)   = tmp_recv2(i,LAP+1:LAP+LAP,3)
                ctr(nx+i,ny+1:ny+LAP)%length(2)   = tmp_recv2(i,LAP+1:LAP+LAP,4)
                ctr(nx+i,ny+1:ny+LAP)%area        = tmp_recv2(i,LAP+1:LAP+LAP,5)
            enddo
        endif
        deallocate(tmp_send1,tmp_send2,tmp_recv1,tmp_recv2)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
        
        !------------------------------------------
        !>J direction
        !------------------------------------------
        allocate(tmp_send1(2*LAP,LAP,Ndata_send),tmp_send2(2*LAP,LAP,Ndata_send))
        allocate(tmp_recv1(2*LAP,LAP,Ndata_send),tmp_recv2(2*LAP,LAP,Ndata_send))
        do j=1,LAP
            tmp_send1(1:LAP,j,1)=ctr(1-LAP:0,j)%x
            tmp_send1(1:LAP,j,2)=ctr(1-LAP:0,j)%y
            tmp_send1(1:LAP,j,3)=ctr(1-LAP:0,j)%length(1)
            tmp_send1(1:LAP,j,4)=ctr(1-LAP:0,j)%length(2)
            tmp_send1(1:LAP,j,5)=ctr(1-LAP:0,j)%area
            tmp_send1(LAP+1:LAP+LAP,j,1)=ctr(nx+1:nx+LAP,j)%x
            tmp_send1(LAP+1:LAP+LAP,j,2)=ctr(nx+1:nx+LAP,j)%y
            tmp_send1(LAP+1:LAP+LAP,j,3)=ctr(nx+1:nx+LAP,j)%length(1)
            tmp_send1(LAP+1:LAP+LAP,j,4)=ctr(nx+1:nx+LAP,j)%length(2)
            tmp_send1(LAP+1:LAP+LAP,j,5)=ctr(nx+1:nx+LAP,j)%area
                
            tmp_send2(1:LAP,j,1)=ctr(1-LAP:0,ny+j-LAP)%x
            tmp_send2(1:LAP,j,2)=ctr(1-LAP:0,ny+j-LAP)%y
            tmp_send2(1:LAP,j,3)=ctr(1-LAP:0,ny+j-LAP)%length(1)
            tmp_send2(1:LAP,j,4)=ctr(1-LAP:0,ny+j-LAP)%length(2)
            tmp_send2(1:LAP,j,5)=ctr(1-LAP:0,ny+j-LAP)%area
            tmp_send2(LAP+1:LAP+LAP,j,1)=ctr(nx+1:nx+LAP,ny+j-LAP)%x
            tmp_send2(LAP+1:LAP+LAP,j,2)=ctr(nx+1:nx+LAP,ny+j-LAP)%y
            tmp_send2(LAP+1:LAP+LAP,j,3)=ctr(nx+1:nx+LAP,ny+j-LAP)%length(1)
            tmp_send2(LAP+1:LAP+LAP,j,4)=ctr(nx+1:nx+LAP,ny+j-LAP)%length(2)
            tmp_send2(LAP+1:LAP+LAP,j,5)=ctr(nx+1:nx+LAP,ny+j-LAP)%area
        enddo

        call MPI_Sendrecv(tmp_send1, LAP*2*LAP*Ndata_send, MPI_KREAL, ID_Y_minus, 9000, &
                          tmp_recv2, LAP*2*LAP*Ndata_send, MPI_KREAL, ID_Y_plus,  9000, MPI_COMM_WORLD, Status, ierr)
        call MPI_Sendrecv(tmp_send2, LAP*2*LAP*Ndata_send, MPI_KREAL, ID_Y_plus,  8000,   &
                          tmp_recv1, LAP*2*LAP*Ndata_send, MPI_KREAL, ID_Y_minus, 8000, MPI_COMM_WORLD, Status, ierr)


        if(ID_Y_minus .ne. MPI_PROC_NULL) then
            do j=1,LAP
                ctr(1-LAP:0,j-LAP)%x          = tmp_recv1(1:LAP,j,1)
                ctr(1-LAP:0,j-LAP)%y          = tmp_recv1(1:LAP,j,2)
                ctr(1-LAP:0,j-LAP)%length(1)  = tmp_recv1(1:LAP,j,3)
                ctr(1-LAP:0,j-LAP)%length(2)  = tmp_recv1(1:LAP,j,4)
                ctr(1-LAP:0,j-LAP)%area       = tmp_recv1(1:LAP,j,5)
                ctr(nx+1:nx+LAP,j-LAP)%x          = tmp_recv1(LAP+1:LAP+LAP,j,1)
                ctr(nx+1:nx+LAP,j-LAP)%y          = tmp_recv1(LAP+1:LAP+LAP,j,2)
                ctr(nx+1:nx+LAP,j-LAP)%length(1)  = tmp_recv1(LAP+1:LAP+LAP,j,3)
                ctr(nx+1:nx+LAP,j-LAP)%length(2)  = tmp_recv1(LAP+1:LAP+LAP,j,4)
                ctr(nx+1:nx+LAP,j-LAP)%area       = tmp_recv1(LAP+1:LAP+LAP,j,5)
            enddo
        endif

        if(ID_Y_plus .ne. MPI_PROC_NULL) then
            do j=1,LAP
                ctr(1-LAP:0,ny+j)%x           = tmp_recv2(1:LAP,j,1)
                ctr(1-LAP:0,ny+j)%y           = tmp_recv2(1:LAP,j,2)
                ctr(1-LAP:0,ny+j)%length(1)   = tmp_recv2(1:LAP,j,3)
                ctr(1-LAP:0,ny+j)%length(2)   = tmp_recv2(1:LAP,j,4)
                ctr(1-LAP:0,ny+j)%area        = tmp_recv2(1:LAP,j,5)
                ctr(nx+1:nx+LAP,ny+j)%x           = tmp_recv2(LAP+1:LAP+LAP,j,1)
                ctr(nx+1:nx+LAP,ny+j)%y           = tmp_recv2(LAP+1:LAP+LAP,j,2)
                ctr(nx+1:nx+LAP,ny+j)%length(1)   = tmp_recv2(LAP+1:LAP+LAP,j,3)
                ctr(nx+1:nx+LAP,ny+j)%length(2)   = tmp_recv2(LAP+1:LAP+LAP,j,4)
                ctr(nx+1:nx+LAP,ny+j)%area        = tmp_recv2(LAP+1:LAP+LAP,j,5)
            enddo
        endif
        deallocate(tmp_send1,tmp_send2,tmp_recv1,tmp_recv2)
        
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    end subroutine MPIUpdate_Conner_Coordinates
    
    
    subroutine MPIUpdate_Conner_conVars()
        implicit none
        integer                     :: i,j,Ndata_send
        real(KREAL),allocatable     :: tmp_send1(:,:,:),tmp_send2(:,:,:),tmp_recv1(:,:,:),tmp_recv2(:,:,:)

        !Sendrecv: ctr%ConVars(1:4) 
        Ndata_send=4
        
        !------------------------------------------
        !>I direction
        !------------------------------------------
        allocate(tmp_send1(LAP,LAP*2,Ndata_send),tmp_send2(LAP,LAP*2,Ndata_send))
        allocate(tmp_recv1(LAP,LAP*2,Ndata_send),tmp_recv2(LAP,LAP*2,Ndata_send))
        
        
        do i=1,LAP
            do j=1-LAP,0
                tmp_send1(i,j+LAP,1:4)=ctr(i,j)%ConVars(1:4)
                tmp_send2(i,j+LAP,1:4)=ctr(nx-LAP+i,j)%ConVars(1:4)
            enddo
            do j=ny+1,ny+LAP
                tmp_send1(i,j-ny+LAP,1:4)=ctr(i,j)%ConVars(1:4)
                tmp_send2(i,j-ny+LAP,1:4)=ctr(nx-LAP+i,j)%ConVars(1:4)
            enddo
        enddo
        

        call MPI_Sendrecv(tmp_send1, LAP*LAP*2*Ndata_send, MPI_KREAL, ID_X_minus, 7000, &
                          tmp_recv2, LAP*LAP*2*Ndata_send, MPI_KREAL, ID_X_plus,  7000, MPI_COMM_WORLD, Status, ierr)
        call MPI_Sendrecv(tmp_send2, LAP*LAP*2*Ndata_send, MPI_KREAL, ID_X_plus,  6000,   &
                          tmp_recv1, LAP*LAP*2*Ndata_send, MPI_KREAL, ID_X_minus, 6000, MPI_COMM_WORLD, Status, ierr)

        if(ID_X_minus .ne. MPI_PROC_NULL) then
            do i=1,LAP
                do j=1-LAP,0
                    ctr(i-LAP,j)%ConVars(1:4) = tmp_recv1(i,j+LAP,1:4)
                enddo
                do j=ny+1,ny+LAP
                    ctr(i-LAP,j)%ConVars(1:4) = tmp_recv1(i,j-ny+LAP,1:4)
                enddo
            enddo
        endif

        if(ID_X_plus .ne. MPI_PROC_NULL) then 
            do i=1,LAP
                do j=1-LAP,0
                    ctr(nx+i,j)%ConVars(1:4) = tmp_recv2(i,j+LAP,1:4)
                enddo
                do j=ny+1,ny+LAP
                    ctr(nx+i,j)%ConVars(1:4) = tmp_recv2(i,j-ny+LAP,1:4)
                enddo
            enddo
        endif
        deallocate(tmp_send1,tmp_send2,tmp_recv1,tmp_recv2)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
        
        !------------------------------------------
        !>J direction
        !------------------------------------------
        allocate(tmp_send1(2*LAP,LAP,Ndata_send),tmp_send2(2*LAP,LAP,Ndata_send))
        allocate(tmp_recv1(2*LAP,LAP,Ndata_send),tmp_recv2(2*LAP,LAP,Ndata_send))
        
        do j=1,LAP
            do i=1-LAP,0
                tmp_send1(i+LAP,j,1:4)=ctr(i,j)%ConVars(1:4)
                tmp_send2(i+LAP,j,1:4)=ctr(i,ny+j-LAP)%ConVars(1:4)
            enddo
            do i=nx+1,nx+LAP
                tmp_send1(i-nx+LAP,j,1:4)=ctr(i,j)%ConVars(1:4)
                tmp_send2(i-nx+LAP,j,1:4)=ctr(i,ny+j-LAP)%ConVars(1:4)
            enddo
        enddo

        call MPI_Sendrecv(tmp_send1, LAP*2*LAP*Ndata_send, MPI_KREAL, ID_Y_minus, 9000, &
                          tmp_recv2, LAP*2*LAP*Ndata_send, MPI_KREAL, ID_Y_plus,  9000, MPI_COMM_WORLD, Status, ierr)
        call MPI_Sendrecv(tmp_send2, LAP*2*LAP*Ndata_send, MPI_KREAL, ID_Y_plus,  8000,   &
                          tmp_recv1, LAP*2*LAP*Ndata_send, MPI_KREAL, ID_Y_minus, 8000, MPI_COMM_WORLD, Status, ierr)


        if(ID_Y_minus .ne. MPI_PROC_NULL) then
            do j=1,LAP
                do i=1-LAP,0
                    ctr(i,j-LAP)%ConVars(1:4) = tmp_recv1(i+LAP,j,1:4)
                enddo
                do i=nx+1,nx+LAP
                    ctr(i,j-LAP)%ConVars(1:4) = tmp_recv1(i-nx+LAP,j,1:4)
                enddo
            enddo
        endif

        if(ID_Y_plus .ne. MPI_PROC_NULL) then
            do j=1,LAP
                do i=1-LAP,0
                    ctr(i,ny+j)%ConVars(1:4) = tmp_recv2(i+LAP,j,1:4)
                enddo
                do i=nx+1,nx+LAP
                    ctr(i,ny+j)%ConVars(1:4) = tmp_recv2(i-nx+LAP,j,1:4)
                enddo
            enddo
        endif
        deallocate(tmp_send1,tmp_send2,tmp_recv1,tmp_recv2)
        
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    end subroutine MPIUpdate_Conner_conVars
    
    
    subroutine MPIUpdate_Conner_VDF()
        implicit none
        integer                                     :: i,j,Ndata_send
        real(KREAL),allocatable,dimension(:,:,:,:)  :: tmp_send1,tmp_send2,tmp_recv1,tmp_recv2

        !Sendrecv: flattened cell-center VDFs h_vec/b_vec
        Ndata_send=2
        
        !------------------------------------------
        !>I direction
        !------------------------------------------
        allocate(tmp_send1(LAP,LAP*2,Nvel,Ndata_send),tmp_send2(LAP,LAP*2,Nvel,Ndata_send))
        allocate(tmp_recv1(LAP,LAP*2,Nvel,Ndata_send),tmp_recv2(LAP,LAP*2,Nvel,Ndata_send))
        do i=1,LAP
            do j=1-LAP,0
                tmp_send1(i,j+LAP,1:Nvel,1)=ctr(i,j)%h_vec(1:Nvel)
                tmp_send1(i,j+LAP,1:Nvel,2)=ctr(i,j)%b_vec(1:Nvel)
                tmp_send2(i,j+LAP,1:Nvel,1)=ctr(nx-LAP+i,j)%h_vec(1:Nvel)
                tmp_send2(i,j+LAP,1:Nvel,2)=ctr(nx-LAP+i,j)%b_vec(1:Nvel)
            enddo
            do j=ny+1,ny+LAP
                tmp_send1(i,j-ny+LAP,1:Nvel,1)=ctr(i,j)%h_vec(1:Nvel)
                tmp_send1(i,j-ny+LAP,1:Nvel,2)=ctr(i,j)%b_vec(1:Nvel)
                tmp_send2(i,j-ny+LAP,1:Nvel,1)=ctr(nx-LAP+i,j)%h_vec(1:Nvel)
                tmp_send2(i,j-ny+LAP,1:Nvel,2)=ctr(nx-LAP+i,j)%b_vec(1:Nvel)
            enddo   
        enddo

        call MPI_Sendrecv(tmp_send1, LAP*LAP*2*Nvel*Ndata_send, MPI_KREAL, ID_X_minus, 7000, &
                          tmp_recv2, LAP*LAP*2*Nvel*Ndata_send, MPI_KREAL, ID_X_plus,  7000, MPI_COMM_WORLD, Status, ierr)
        call MPI_Sendrecv(tmp_send2, LAP*LAP*2*Nvel*Ndata_send, MPI_KREAL, ID_X_plus,  6000,   &
                          tmp_recv1, LAP*LAP*2*Nvel*Ndata_send, MPI_KREAL, ID_X_minus, 6000, MPI_COMM_WORLD, Status, ierr)

        if(ID_X_minus .ne. MPI_PROC_NULL) then
            do i=1,LAP
                do j=1-LAP,0
                    ctr(i-LAP,j)%h_vec(1:Nvel) = tmp_recv1(i,j+LAP,1:Nvel,1)
                    ctr(i-LAP,j)%b_vec(1:Nvel) = tmp_recv1(i,j+LAP,1:Nvel,2)
                enddo
                do j=ny+1,ny+LAP
                    ctr(i-LAP,j)%h_vec(1:Nvel) = tmp_recv1(i,j-ny+LAP,1:Nvel,1)
                    ctr(i-LAP,j)%b_vec(1:Nvel) = tmp_recv1(i,j-ny+LAP,1:Nvel,2)
                enddo
            enddo
        endif

        if(ID_X_plus .ne. MPI_PROC_NULL) then
            do i=1,LAP
                do j=1-LAP,0
                    ctr(nx+i,j)%h_vec(1:Nvel) = tmp_recv2(i,j+LAP,1:Nvel,1)
                    ctr(nx+i,j)%b_vec(1:Nvel) = tmp_recv2(i,j+LAP,1:Nvel,2)
                enddo
                do j=ny+1,ny+LAP
                    ctr(nx+i,j)%h_vec(1:Nvel) = tmp_recv2(i,j-ny+LAP,1:Nvel,1)
                    ctr(nx+i,j)%b_vec(1:Nvel) = tmp_recv2(i,j-ny+LAP,1:Nvel,2)
                enddo
            enddo
        endif
        deallocate(tmp_send1,tmp_send2,tmp_recv1,tmp_recv2)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
        
        !------------------------------------------
        !>J direction
        !------------------------------------------
        allocate(tmp_send1(2*LAP,LAP,Nvel,Ndata_send),tmp_send2(2*LAP,LAP,Nvel,Ndata_send))
        allocate(tmp_recv1(2*LAP,LAP,Nvel,Ndata_send),tmp_recv2(2*LAP,LAP,Nvel,Ndata_send))
        do j=1,LAP
            do i=1-LAP,0
                tmp_send1(i+LAP,j,1:Nvel,1)=ctr(i,j)%h_vec(1:Nvel)
                tmp_send1(i+LAP,j,1:Nvel,2)=ctr(i,j)%b_vec(1:Nvel)
                tmp_send2(i+LAP,j,1:Nvel,1)=ctr(i,ny+j-LAP)%h_vec(1:Nvel)
                tmp_send2(i+LAP,j,1:Nvel,2)=ctr(i,ny+j-LAP)%b_vec(1:Nvel)
            enddo
            do i=nx+1,nx+LAP
                tmp_send1(i-nx+LAP,j,1:Nvel,1)=ctr(i,j)%h_vec(1:Nvel)
                tmp_send1(i-nx+LAP,j,1:Nvel,2)=ctr(i,j)%b_vec(1:Nvel)
                tmp_send2(i-nx+LAP,j,1:Nvel,1)=ctr(i,ny+j-LAP)%h_vec(1:Nvel)
                tmp_send2(i-nx+LAP,j,1:Nvel,2)=ctr(i,ny+j-LAP)%b_vec(1:Nvel)
            enddo
        enddo

        call MPI_Sendrecv(tmp_send1, LAP*2*LAP*Nvel*Ndata_send, MPI_KREAL, ID_Y_minus, 9000, &
                          tmp_recv2, LAP*2*LAP*Nvel*Ndata_send, MPI_KREAL, ID_Y_plus,  9000, MPI_COMM_WORLD, Status, ierr)
        call MPI_Sendrecv(tmp_send2, LAP*2*LAP*Nvel*Ndata_send, MPI_KREAL, ID_Y_plus,  8000,   &
                          tmp_recv1, LAP*2*LAP*Nvel*Ndata_send, MPI_KREAL, ID_Y_minus, 8000, MPI_COMM_WORLD, Status, ierr)


        if(ID_Y_minus .ne. MPI_PROC_NULL) then
            do j=1,LAP
                do i=1-LAP,0
                    ctr(i,j-LAP)%h_vec(1:Nvel) = tmp_recv1(i+LAP,j,1:Nvel,1)
                    ctr(i,j-LAP)%b_vec(1:Nvel) = tmp_recv1(i+LAP,j,1:Nvel,2)
                enddo
                do i=nx+1,nx+LAP
                    ctr(i,j-LAP)%h_vec(1:Nvel) = tmp_recv1(i-nx+LAP,j,1:Nvel,1)
                    ctr(i,j-LAP)%b_vec(1:Nvel) = tmp_recv1(i-nx+LAP,j,1:Nvel,2)
                enddo
            enddo
        endif

        if(ID_Y_plus .ne. MPI_PROC_NULL) then
            do j=1,LAP
                do i=1-LAP,0
                    ctr(i,ny+j)%h_vec(1:Nvel) = tmp_recv2(i+LAP,j,1:Nvel,1)
                    ctr(i,ny+j)%b_vec(1:Nvel) = tmp_recv2(i+LAP,j,1:Nvel,2)
                enddo
                do i=nx+1,nx+LAP
                    ctr(i,ny+j)%h_vec(1:Nvel) = tmp_recv2(i-nx+LAP,j,1:Nvel,1)
                    ctr(i,ny+j)%b_vec(1:Nvel) = tmp_recv2(i-nx+LAP,j,1:Nvel,2)
                enddo
            enddo
        endif
        deallocate(tmp_send1,tmp_send2,tmp_recv1,tmp_recv2)
        
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
     end subroutine MPIUpdate_Conner_VDF
     
    
    !------------------------------------------------------------------------------------------------
    !>The following subroutines are more general, but only support sendrecv one variables at one time
    !------------------------------------------------------------------------------------------------
    !Send and recv mesg: updata overlaped macro-variables in x direction (one variables at one time)
    subroutine MPIUpdata_MacroVars_Xdirc(f)
        implicit none
        real(KREAL),intent(inout)   :: f(1-LAP:nx+LAP,1-LAP:ny+LAP)
        integer                     :: i,j,k1
        real(KREAL)                 :: tmp_send1(LAP*ny),tmp_send2(LAP*ny),tmp_recv1(LAP*ny),tmp_recv2(LAP*ny)

        do j=1,ny
            do i=1,LAP
                k1=(j-1)*LAP+i
                tmp_send1(k1)=f(i,j)
                tmp_send2(k1)=f(nx-LAP+i,j)
            enddo
        enddo

        call MPI_Sendrecv(tmp_send1, LAP*ny, MPI_KREAL, ID_X_minus, 9000, &
                          tmp_recv2, LAP*ny, MPI_KREAL, ID_X_plus,  9000, MPI_COMM_WORLD, Status, ierr)
        call MPI_Sendrecv(tmp_send2, LAP*ny, MPI_KREAL, ID_X_plus,  8000,   &
                          tmp_recv1, LAP*ny, MPI_KREAL, ID_X_minus, 8000, MPI_COMM_WORLD, Status, ierr)

        if(ID_X_minus .ne. MPI_PROC_NULL) then
            do j=1,ny
                do i=1,LAP
                    k1=(j-1)*LAP+i
                    f(i-LAP,j)=tmp_recv1(k1)
                enddo
            enddo
        endif

        if(ID_X_plus .ne. MPI_PROC_NULL) then
            do j=1,ny
                do i=1,LAP
                    k1=(j-1)*LAP+i
                    f(nx+i,j)=tmp_recv2(k1)
                enddo
            enddo
        endif

    end subroutine MPIUpdata_MacroVars_Xdirc

    
    !Send and recv mesg: updata overlaped macro-variables in y direction (one variables at one time)
    subroutine MPIUpdata_MacroVars_Ydirc(f)
        implicit none
        real(KREAL),intent(inout)   :: f(1-LAP:nx+LAP,1-LAP:ny+LAP)
        integer                     :: i,j,k1
        real(KREAL)                 :: tmp_send1(LAP*nx),tmp_send2(LAP*nx),tmp_recv1(LAP*nx),tmp_recv2(LAP*nx)

        do j=1,LAP
            do i=1,nx
                k1=(j-1)*nx+i
                tmp_send1(k1)=f(i,j)
                tmp_send2(k1)=f(i,ny+j-LAP)
            enddo
        enddo

        call MPI_Sendrecv(tmp_send1, LAP*nx, MPI_KREAL, ID_Y_minus, 9000, &
                          tmp_recv2, LAP*nx, MPI_KREAL, ID_Y_plus,  9000, MPI_COMM_WORLD, Status, ierr)
        call MPI_Sendrecv(tmp_send2, LAP*nx, MPI_KREAL, ID_Y_plus,  8000,   &
                          tmp_recv1, LAP*nx, MPI_KREAL, ID_Y_minus, 8000, MPI_COMM_WORLD, Status, ierr)


        if(ID_Y_minus .ne. MPI_PROC_NULL) then
            do j=1,LAP
                do i=1,nx
                    k1=(j-1)*nx+i
                    f(i,j-LAP)=tmp_recv1(k1)
                enddo
            enddo
        endif

        if(ID_Y_plus .ne. MPI_PROC_NULL) then
            do j=1,LAP
                do i=1,nx
                    k1=(j-1)*nx+i
                    f(i,ny+j)=tmp_recv2(k1)
                enddo
            enddo
        endif

    end subroutine MPIUpdata_MacroVars_Ydirc

    !Exchange one face-neighbour halo layer for a packed correction field.
    subroutine MPIExchangeField3D(field,ncomp,tag_base)
        integer(KINT), intent(in) :: ncomp,tag_base
        real(KREAL), intent(inout) :: field(:,1-LAP:,1-LAP:)
        real(KREAL), allocatable :: sxm(:),sxp(:),rxm(:),rxp(:)
        real(KREAL), allocatable :: sym(:),syp(:),rym(:),ryp(:)
        integer(KINT) :: i,j,m,p

        allocate(sxm(ncomp*ny),sxp(ncomp*ny),rxm(ncomp*ny),rxp(ncomp*ny))
        p=0
        do j=1,ny
            do m=1,ncomp
                p=p+1; sxm(p)=field(m,1,j); sxp(p)=field(m,nx,j)
            enddo
        enddo
        call MPI_Sendrecv(sxm,size(sxm),MPI_KREAL,ID_X_minus,tag_base, &
                          rxp,size(rxp),MPI_KREAL,ID_X_plus,tag_base, &
                          MPI_COMM_WORLD,Status,ierr)
        call MPI_Sendrecv(sxp,size(sxp),MPI_KREAL,ID_X_plus,tag_base+1, &
                          rxm,size(rxm),MPI_KREAL,ID_X_minus,tag_base+1, &
                          MPI_COMM_WORLD,Status,ierr)
        if (ID_X_minus/=MPI_PROC_NULL) then
            p=0
            do j=1,ny
                do m=1,ncomp
                    p=p+1; field(m,0,j)=rxm(p)
                enddo
            enddo
        endif
        if (ID_X_plus/=MPI_PROC_NULL) then
            p=0
            do j=1,ny
                do m=1,ncomp
                    p=p+1; field(m,nx+1,j)=rxp(p)
                enddo
            enddo
        endif
        deallocate(sxm,sxp,rxm,rxp)

        allocate(sym(ncomp*nx),syp(ncomp*nx),rym(ncomp*nx),ryp(ncomp*nx))
        p=0
        do i=1,nx
            do m=1,ncomp
                p=p+1; sym(p)=field(m,i,1); syp(p)=field(m,i,ny)
            enddo
        enddo
        call MPI_Sendrecv(sym,size(sym),MPI_KREAL,ID_Y_minus,tag_base+2, &
                          ryp,size(ryp),MPI_KREAL,ID_Y_plus,tag_base+2, &
                          MPI_COMM_WORLD,Status,ierr)
        call MPI_Sendrecv(syp,size(syp),MPI_KREAL,ID_Y_plus,tag_base+3, &
                          rym,size(rym),MPI_KREAL,ID_Y_minus,tag_base+3, &
                          MPI_COMM_WORLD,Status,ierr)
        if (ID_Y_minus/=MPI_PROC_NULL) then
            p=0
            do i=1,nx
                do m=1,ncomp
                    p=p+1; field(m,i,0)=rym(p)
                enddo
            enddo
        endif
        if (ID_Y_plus/=MPI_PROC_NULL) then
            p=0
            do i=1,nx
                do m=1,ncomp
                    p=p+1; field(m,i,ny+1)=ryp(p)
                enddo
            enddo
        endif
        deallocate(sym,syp,rym,ryp)
    end subroutine MPIExchangeField3D

    subroutine MPIExchangeDeltaMacro(field)
        use mod_ImplicitState, only: implicit_exchange_seconds
        real(KREAL), intent(inout) :: field(:,1-LAP:,1-LAP:)
        real(KREAL) :: start_time

        start_time=MPI_Wtime()
        call MPIExchangeField3D(field,4,17000)
        implicit_exchange_seconds=implicit_exchange_seconds+MPI_Wtime()-start_time
    end subroutine MPIExchangeDeltaMacro

    subroutine MPIExchangeDeltaMicro(field_h,field_b)
        use mod_ImplicitState, only: implicit_exchange_seconds
        real(KREAL), intent(inout) :: field_h(:,1-LAP:,1-LAP:),field_b(:,1-LAP:,1-LAP:)
        real(KREAL) :: start_time

        start_time=MPI_Wtime()
        call MPIExchangeField3D(field_h,Nvel,18000)
        call MPIExchangeField3D(field_b,Nvel,18100)
        implicit_exchange_seconds=implicit_exchange_seconds+MPI_Wtime()-start_time
    end subroutine MPIExchangeDeltaMicro

end module mpi_wrapper
