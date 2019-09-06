subroutine passe(l, edge, parent, edge_length, n_species, traits, n_splits1, n_splits2, &
                   mu0, mu1, q01, q10, l0, l1, n_pars, diag)

  implicit none
  integer, intent(IN) :: n_species, n_splits1, n_splits2, n_pars
  integer :: i, m, n, n_splits
  integer, dimension(2*n_species - 2), intent(IN) :: edge
  integer, dimension(2*n_species - 1), intent(IN) :: parent
  integer, dimension(2*n_species - 2) :: done
  
  double precision, dimension(n_pars), intent(IN) :: mu0, mu1, q01, q10, l0, l1
  double precision, dimension(n_pars) :: l, e0c, e1c, dn0c, dn1c, a, b, aa, bb, &
    d1, d2, dk, c1, c2, cc, k1, k2, ek1, ek2, x, cc1, cc2
  double precision ::  tn, tc, g, tnmax
  double precision, dimension(2*n_species - 2), intent(IN) :: edge_length
  double precision, dimension(2*n_species - 2, n_pars) :: e0, e1, dn0, dn1
  double precision, dimension(2*n_species - 1, n_pars) :: lq
  double precision, dimension(n_species), intent(IN) :: traits
  
  logical :: diag
  
  ! diagnostic file
  if (diag) then
    open(10, file = 'out.txt')
  endif
  
  ! nothing done yet
  done = 0
  
  ! initialise variables
  e0 = 0.
  e1 = 0.
  lq = 0.

  ! precalculations
  !!$acc kernels
  a = mu0 + q01 + l0
  b = mu1 + q10 + l1
  !!$acc end kernels
  
  if (diag) then
    write(10, *) "Populate known traits"
    flush(10)
  endif
  
  ! populate species with known traits with probabilities
  do i = 1,(2*n_species - 2)
    if (edge(i) <= n_species) then
      !!$acc kernels
      dn0(i,:) = 1. - traits(edge(i))
      dn1(i,:) = traits(edge(i))
      !!$acc end kernels
      done(i) = 1
    endif
  enddo
  
  if (diag) then
    write(10, *) "Populate all other species"
    flush(10)
  endif
  
  ! populate all other species by integration
  do while (any(done < 2))
    if (diag) then
      write(10, *) "Scan for populated children"
      flush(10)
    endif
    
    do i = 1, (n_species - 1)
    ! search for parents with populated children
    if (done(2*i - 1) .eq. 1 .and. done(2*i) .eq. 1) then
      done(2*i - 1) = 2
      done(2*i) = 2

      if (diag) then
        write(10, *) "Integrate children", i
        flush(10)
      endif
      ! estimate updated parameters at parent
      ! starting from each child separately
      do m = 0, 1
        if (diag) then
          write(10, *) "Integrate child", m
          flush(10)
        endif
        
        e0c = e0(2*i - 1 + m, :)
        e1c = e1(2*i - 1 + m, :)
        dn0c = dn0(2*i - 1 + m, :)
        dn1c = dn1(2*i - 1 + m, :)
        tc = 0.
        
        if (edge(2*i - 1 + m) <= n_species) then
          n_splits = n_splits1
        else
          n_splits = n_splits2
        endif
        
        ! integrate in n_splits steps
        do n = 1, n_splits
          if (diag) then
            write(10, *) "Integrate step", n
            write(10, *) "Start", e0c, e1c, dn0c, dn1c
            flush(10)
          endif
          !!$acc kernels
          d1 = mu0 - a*e0c + q01*e1c + l0*e0c*e0c
          aa = a - 2.*l0*e0c
          d2 = mu1 - b*e1c + q10*e0c + l1*e1c*e1c
          bb = b - 2.*l1*e1c
          !!$acc end kernels
          !!$acc kernels
          cc = aa + bb
          !!$acc end kernels
          !!$acc kernels
          dk = sqrt(cc*cc - 4.*(aa*bb - q01*q10))
          !!$acc end kernels
          !!$acc kernels
          k1 = (cc + dk)/2.
          k2 = (cc - dk)/2.
          !!$acc end kernels

          tnmax = edge_length(2*i - 1 + m) - tc
          
          if (n < n_splits) then
            tn = edge_length(2*i - 1 + m)/(n_splits + 0.)
            !tn = -log(1. - 1./(n_splits - n + 1.))/max(k1,k2)
            if (tn > tnmax) then
              tn = tnmax
            endif
          else
            tn = tnmax
          endif
          
          if (diag) then
            write(10, *) "Time", tn
            flush(10)
          endif
          !!$acc kernels
          ek1 = 1. - exp(-tn*k1)
          ek2 = 1. - exp(-tn*k2)
          c1 = -((k2 - aa)*d1 + q01*d2)/dk/k1
          c2 = ((k1 - aa)*d1 + q01*d2)/dk/k2
          e0c = e0c + c1*ek1 + c2*ek2
          e1c = e1c + c1*ek1*(aa-k1)/q01 + c2*ek2*(aa-k2)/q01

          d1 = -aa*dn0c + q01*dn1c
          d2 = -bb*dn1c + q10*dn0c
          cc1 = -((k2 - aa)*d1 + q01*d2)/dk/k1
          cc2 = ((k1 - aa)*d1 + q01*d2)/dk/k2
          dn0c = dn0c + cc1*ek1 + cc2*ek2
          dn1c = dn1c + cc1*ek1*(aa-k1)/q01 + cc2*ek2*(aa-k2)/q01

          if (diag) then
            write(10,*) "k", k1, k2
            write(10,*) "c", c1, c2
            write(10,*) "cc", cc1, cc2
            write(10,*) "mod1", (aa-k1)/q01
            write(10,*) "mod2", (aa-k2)/q01
            flush(10)
          endif

          !!$acc end kernels
          tc = tc + tn
          if (diag) then
            write(10, *) "End", e0c, e1c, dn0c, dn1c
            flush(10)
          endif
          ! if reached the end then skip out
          if (tc >= edge_length(2*i - 1 + m)) then
            exit
          endif
      enddo ! end current split
    
    !!$acc kernels
      e0(2*i - 1 + m, :) = e0c
      e1(2*i - 1 + m, :) = e1c
      x = dn0c + dn1c
      dn0(2*i - 1 + m, :) = dn0c/x
      dn1(2*i - 1 + m, :) = dn1c/x
      lq(edge(2*i - 1 + m), :) = lq(edge(2*i - 1 + m), :) + log(x)
      !!$acc end kernels
    enddo ! end current child
    !!$acc kernels
    dn0c = dn0(2*i - 1, :)*dn0(2*i, :)*l0
    dn1c = dn1(2*i - 1, :)*dn1(2*i, :)*l1
    x = dn0c + dn1c
    dn0c = dn0c / x
    dn1c = dn1c / x
    !!$acc end kernels
    if (i .ne. 1) then
      lq(n_species + i, :) = lq(n_species + i, :) + log(x)
    
      ! populate parent based on both child results
      !!$acc kernels
      e0(parent(i), :) = max(e0(2*i - 1, :), e0(2*i, :))
      e1(parent(i), :) = max(e1(2*i - 1, :), e1(2*i, :))
      dn0(parent(i), :) = dn0c
      dn1(parent(i), :) = dn1c
      done(parent(i)) = 1
      !!$acc end kernels
    endif
    
    if (diag) then
      write(10, *) "Populate parent"
      flush(10)
    endif
    
    !endif ! if false
    
    ! unnecessary code left just in case
    !where (edge(:,2) .eq. (n_species + i))
    !end where
    endif
    enddo ! end current branch
  enddo
  
  
  ! calculate values at root
  !!$acc kernels
  dn0c = dn0(1, :)*dn0(2, :)*l0
  dn1c = dn1(1, :)*dn1(2, :)*l1
  !!$acc end kernels
  
  ! calculate final likelihood 
  if (diag) then
    write(10, *) "Final calculation"
    flush(10)
  endif
  
  ! weighting of root states based on
  ! equilibrium state frequencies
  ! (Maddinson 2007, Appendix 2)
  !  g = l0 + mu1 - l1 - mu0
  !  if (g .eq. 0.) then
  !    x = q10 / (q01 + q10)
  !  else
    !    b = g - q01 - q10
  !    x = (b + sqrt(b*b + 4.*g*q10))/2./g
  !  endif
  
  ! FitzJohn et al. 2009
  !!$acc kernels
  x = dn0c / (dn0c + dn1c)
  l = dn0c*x + dn1c*(1.-x)
  
  if (diag) then
    write(10, *) "x", x
    write(10, *) "l", l
    flush(10)
  endif
  
  l = sum(lq, 1) + log(l) - log(x*l0*(1.-e0(1, :))*(1.-e0(1, :)) + &
        (1.-x)*l1*(1.-e1(1, :))*(1.-e1(1, :)))
  !!$acc end kernels
  if (diag) then
    write(10, *) "l2", l
    close(10)
  endif

end subroutine


subroutine fast_bisse(x0, n_max, n_species, edge, edge_length, edge_state, &
                    node_state, tip_state, mu0, mu1, q01, q10, l0, l1, diag)

  implicit none
  integer :: n0, n1, row, rows, species, curr_parent, curr_sibling, &
             event_state, event_type, idx, n, i, root_state
  integer, intent(IN) :: n_species, n_max
  integer, dimension(2*n_species - 1) :: tipnode
  integer, dimension(2*n_species - 2) :: edge_state, ei
  integer, dimension(2*n_species - 2, 2) :: edge
  integer, dimension(n_species) :: tip_state
  integer, dimension(n_species - 1) :: node_state
  
  integer, dimension(2*n_max - 2, 2) :: edge_tmp
  integer, dimension(2*n_max - 2) :: edge_state_tmp, row0, row1, parent, sibling, labels
  
  double precision :: x0, rand, t, t_event, tmp1, tmp2, tmp
  double precision, intent(IN) :: mu0, mu1, q01, q10, l0, l1
  double precision, dimension(2*n_species - 2) :: edge_length
  double precision, dimension(2*n_max - 2) :: edge_length_tmp
  
  logical :: diag, root
  logical, dimension(2*n_species - 2) :: istip
  logical, dimension(2*n_max - 2) :: isedge, istip_tmp
  
  ! diagnostic file
  if (diag) then
    open(10, file = 'outfast.txt')
  endif

  ! set random seed (note that this currently cannot be user set)  
  call init_random_seed() 
  n = 0

  if (diag) then
    write(10, *) 'Start'
    flush(10)
  endif
  
  call random_number(rand) ! generate uniform random (0,1)
  if ((x0 .ne. 0D0) .and. (x0 .ne. 1D0)) then
    if (rand < x0) then
      root_state = 0
    else
      root_state = 1
    endif
  else
    root_state = int(x0)
  endif
  
  n0 = 0
  n1 = 0

  if (diag) then
    write(10, *) 'Start loop'
    flush(10)
  endif
  
  do while ((n0 + n1) <= n_species .and. n < n_max)
    n = n + 1
    do while ((n0 + n1) < 2)
      isedge = .false.
      t = 0.
      if (root_state == 0) then
        n0 = 2
        n1 = 0
        row0(1:2) = (/ 1, 2 /)
      else
        n0 = 0
        n1 = 2
        row1(1:2) = (/ 1, 2 /)
      endif
    
      edge_tmp(1,:) = (/ 1, 2 /)
      edge_tmp(2,:) = (/ 1, 3 /)
      parent(1:2) = (/ 0, 0 /)
      sibling(1:2) = (/ 2, 1 /)
      edge_state_tmp(1:2) = root_state
      edge_length_tmp = 0.
      rows = 2
      isedge(1:2) = .true.
      species = 4
      if (diag) then
        write(10, *) 'Initialisation done'
        flush(10)
      endif
    enddo
    
    if (diag) then
      write(10, *) 'Species', n0, n1
      flush(10)
    endif
    tmp1 = (l0+mu0+q01)*n0 
    tmp2 = (l1+mu1+q10)*n1
    tmp = tmp1 + tmp2
    call random_number(rand) ! generate uniform random (0,1)
    t_event = -log(rand)/tmp ! time till next event
    
    if (diag) then
      write(10, *) 'Event time', t_event
      flush(10)
    endif
    
    if (n0 > 0) then
      edge_length_tmp(row0(1:n0)) = edge_length_tmp(row0(1:n0)) + t_event
    endif
    if (n1 > 0) then
      edge_length_tmp(row1(1:n1)) = edge_length_tmp(row1(1:n1)) + t_event
    endif
    
    ! if we've reached max species, don't actually work out
    ! the next event, just stop here
    if ((n0 + n1) == n_species) then
      exit
    endif
    
    call random_number(rand) ! generate uniform random (0,1)
    
    ! work out which state has an event, then which species,
    ! then which type of event, and update number of species for each state
    
    if (rand < (tmp1 / tmp) .and. n0 > 0) then
      event_state = 0
      call random_number(rand) ! generate uniform random (0,1)  
      idx = ceiling(rand*n0)
      row = row0(idx)
      call random_number(rand) ! generate uniform random (0,1)  
      tmp = l0 + mu0 + q01
      if (rand < (l0/tmp)) then
        event_type = 0
        n0 = n0 + 1
      else if (rand < ((l0 + mu0)/tmp)) then
        event_type = 1
        n0 = n0 - 1
      else
        event_type = 2
        n0 = n0 - 1
        n1 = n1 + 1
      endif
    else if (rand >= (tmp1 / tmp) .and. n1 > 0) then
      event_state = 1
      call random_number(rand) ! generate uniform random (0,1)
      idx = ceiling(rand*n1)
      row = row1(idx)
      call random_number(rand) ! generate uniform random (0,1)  
      tmp = l1 + mu1 + q10
      if (rand < (l1/tmp)) then
        event_type = 0      
        n1 = n1 + 1
      else if (rand < ((l1 + mu1)/tmp)) then
        event_type = 1      
        n1 = n1 - 1
      else
        event_type = 2
        n0 = n0 + 1
        n1 = n1 - 1
      endif
    endif

    if ((n0 + n1) < 2) then
      cycle
    endif

    if (diag) then
      write(10, *) 'State', event_state, 'Type', event_type
      write(10,*) 'Row', row
      flush(10)
    endif
    
    ! change edge data, parent/sibling data, rows for state 0 and 1,
    ! states, edge lengths, total # rows and # species
    
    if (event_type == 0) then ! speciation
      edge_tmp(rows + 1, :) = (/ edge_tmp(row, 2), species /)
      edge_tmp(rows + 2, :) = (/ edge_tmp(row, 2), species + 1 /)
      edge_length_tmp((rows+1):(rows+2)) = 0.
      parent((rows + 1):(rows + 2)) = row
      sibling((rows + 1):(rows + 2)) = (/ rows + 2, rows + 1 /)
      isedge((rows + 1):(rows + 2)) = .true.
      if (event_state == 0) then
        row0(idx) = rows + 1
        row0(n0) = rows + 2
        edge_state_tmp((rows + 1):(rows + 2)) = 0
      else
        row1(idx) = rows + 1
        row1(n1) = rows + 2
        edge_state_tmp((rows + 1):(rows + 2)) = 1
      endif
      rows = rows + 2
      species = species + 2
    else if (event_type == 1) then ! extinction
      ! stomp edge with extinct species
      isedge(row) = .false.
      curr_parent = parent(row)
      curr_sibling = sibling(row)
      
      root = .false.
      if (curr_parent == 0) then
        root = .true.
        curr_parent = curr_sibling
        if (diag) then
          write(10,*) 'removing root edge'
          flush(10)
        endif
      else
        ! inherit new parents and siblings
        parent(curr_sibling) = parent(curr_parent)
        sibling(curr_sibling) = sibling(curr_parent)
        sibling(sibling(curr_parent)) = curr_sibling
      endif  

      if (diag) then
        write(10,*) 'stomped 1'
        flush(10)
      endif
      
      if (event_state == 0) then
        row0(idx:n0) = row0((idx+1):(n0+1))
      else
        row1(idx:n1) = row1((idx+1):(n1+1))
      endif

      if (root) then
        where (edge_tmp(1:rows, 1) == edge_tmp(curr_sibling, 2))
          parent(1:rows) = 0
        end where
      endif

      if (diag) then
        write(10,*) 'updated row stuff 1'
        flush(10)
      endif
      
      if (.not. root) then
        edge_tmp(curr_sibling, 1) = edge_tmp(curr_parent, 1)
        edge_length_tmp(curr_sibling) = edge_length_tmp(curr_parent) + & 
                                        edge_length_tmp(curr_sibling)
        if (diag) then
          write(10,*) 'updated sibling stuff'
          flush(10)
        endif
      endif
      
      ! stomp edge to original parent
      isedge(curr_parent) = .false.

      if (diag) then
        write(10,*) 'stomped 2'
        flush(10)
      endif

    else ! transfer
      edge_state_tmp(row) = 1 - edge_state_tmp(row)
      if (event_state == 0) then
        row0(idx:n0) = row0((idx+1):(n0+1))
        row1(n1) = row
      else
        row1(idx:n1) = row1((idx+1):(n1+1))
        row0(n0) = row
      endif
    
    endif
    
    if (diag) then
      write(10,'(80I5)') edge_tmp(1:rows, 1)
      write(10,'(80I5)') edge_tmp(1:rows, 2)
      write(10,'(80L5)') isedge(1:rows)
      write(10,'(80F7.3)') edge_length_tmp(1:rows)
      write(10,*) 'row0', row0(1:n0)
      write(10,*) 'row1', row1(1:n1)
      write(10,*) 'parent'
      write(10,'(80I5)') parent(1:rows)
      write(10,*) 'sibling'
      write(10,'(80I5)') sibling(1:rows)
      write(10,*) count(isedge)
      flush(10)
    endif
    
  enddo
  
  edge(:,1) = pack(edge_tmp(:, 1), isedge)
  edge(:,2) = pack(edge_tmp(:, 2), isedge)
  edge_length = pack(edge_length_tmp, isedge)
  
  istip_tmp = .false.
  if (n0 > 0) then
    istip_tmp(row0(1:n0)) = .true.
  endif
  if (n1 > 0) then
    istip_tmp(row1(1:n1)) = .true.
  endif
  istip = pack(istip_tmp, isedge)
  
  species = 0
  n = n_species + 1
  labels = n_species + 1
  do i = 1, (2*n_species - 2)
    if (istip(i)) then
      species = species + 1
      labels(edge(i,2)) = species
      edge(i,2) = species
    else
      n = n + 1
      labels(edge(i,2)) = n
      edge(i,2) = n
    endif
  enddo

  ei = labels(edge(:,1))
!!!$acc kernels
  do i = 1, (2*n_species - 2)
    edge(i,1) = ei(i) !labels(edge(i,1))
  enddo
 !!!$acc end kernels

  tipnode = 0
  edge_state = pack(edge_state_tmp, isedge)
  tipnode(edge(:,2)) = edge_state
  tip_state = tipnode(1:n_species)
  node_state = tipnode((n_species+1):(2*n_species - 1))
  node_state(1) = root_state
  
  if (diag) then
    write(10,*) 'End'
    write(10,*) x0
    write(10,*) n_max
    write(10,*) n_species
    write(10,*) edge(:,1)
    write(10,*) edge(:,2)
    write(10,*) edge_length
    write(10,*) edge_state
    write(10,*) node_state
    write(10,*) tip_state
    write(10,*) mu0
    write(10,*) mu1
    write(10,*) q01
    write(10,*) q10
    write(10,*) l0
    write(10,*) l1
    write(10,*) diag
    close(10)
  endif

end subroutine


! uses GCC GNU example routine for initialising random seed
! note that this is currently not settable by user
! https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html

subroutine init_random_seed()
  use iso_fortran_env, only: int64
  implicit none
  integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt(8), pid
  integer(int64) :: t

  call random_seed(size = n)
  allocate(seed(n))
  ! First try if the OS provides a random number generator
  open(newunit=un, file="/dev/urandom", access="stream", &
       form="unformatted", action="read", status="old", iostat=istat)
  if (istat == 0) then
     read(un) seed
     close(un)
  else
     ! Fallback to XOR:ing the current time and pid. The PID is
     ! useful in case one launches multiple instances of the same
     ! program in parallel.
     call system_clock(t)
     if (t == 0) then
        call date_and_time(values=dt)
        t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
             + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
             + dt(3) * 24_int64 * 60 * 60 * 1000 &
             + dt(5) * 60 * 60 * 1000 &
             + dt(6) * 60 * 1000 + dt(7) * 1000 &
             + dt(8)
     end if
     pid = getpid()
     t = ieor(t, int(pid, kind(t)))
     do i = 1, n
        seed(i) = lcg(t)
     end do
  end if
  call random_seed(put=seed)
contains
  ! This simple PRNG might not be good enough for double precision work, but is
  ! sufficient for seeding a better PRNG.
  function lcg(s)
    integer :: lcg
    integer(int64) :: s
    if (s == 0) then
       s = 104729
    else
       s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))
  end function lcg
end subroutine init_random_seed