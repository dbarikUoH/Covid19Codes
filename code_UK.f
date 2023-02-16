c	***************************************************************
c	HDC model of Covid-19 for UK
c	Hafijur Rahaman & Debashis Barik, SoC, University of Hyderabad
c	Date: 14.02.2023
c	***************************************************************

        implicit real *8(a-h,o-z)
        parameter(ips=1000)	    ! size of lattice in X
        parameter(jps=1000)	    ! size of lattice in Y
        dimension ns(ips,jps)	    ! ns(i,j): Susceptible
        dimension nsc(ips,jps)	    ! nsc(i,j): Susceptible-Child 
        dimension nsy(ips,jps)	    ! nsy(i,j): Susceptible-Young 
        dimension nso(ips,jps)	    ! nso(i,j): Susceptible-Old 
        dimension ne(ips,jps)	    ! ne(i,j):  Exposed
        dimension nec(ips,jps)	    ! nec(i,j): Exposed-Child
        dimension ney(ips,jps)	    ! ney(i,j): Exposed-Young
        dimension neo(ips,jps)	    ! neo(i,j): Exposed-Old
        dimension ni(ips,jps)	    ! ni(i,j):  Infectious
        dimension nic(ips,jps)	    ! nic(i,j): Infectious-Child
        dimension niy(ips,jps)	    ! niy(i,j): Infectious-Young
        dimension nio(ips,jps)	    ! nio(i,j): Infectios-Old
        dimension nr(ips,jps)	    ! nr(i,j):  Recovered
        dimension nrc(ips,jps)	    ! nrc(i,j): Recovered-Child
        dimension nry(ips,jps)	    ! nry(i,j): Recovered-Young
        dimension nro(ips,jps)	    ! nro(i,j): Recovered-Old
        dimension nrr(ips,jps)      ! Exposed to Susceptible  
        dimension dth(ips,jps)	    ! dth(i,j):  Dead 
        dimension dthc(ips,jps)	    ! dthc(i,j): Dead-Child
        dimension dthy(ips,jps)	    ! dthy(i,j): Dead-Young
        dimension dtho(ips,jps)	    ! dth(i,j):  Dead-Old
        dimension vir(ips,jps)      ! vir(i,j): virus
        dimension vir0(ips,jps)     ! vir0(i,j): virus@t=0
        dimension itm(ips,jps)	    ! Time since Exposed
        dimension itmc(ips,jps)     ! Time since Exposed-Child
        dimension itmy(ips,jps)	    ! Time since Exposed-Young
        dimension itmo(ips,jps)	    ! Time since Exposed-Old
        dimension nip(ips,jps)      ! New Infectious per day
        dimension nipc(ips,jps)     ! New Child-Infectious per day
        dimension nipy(ips,jps)     ! New Young-Infectious per day
        dimension nipo(ips,jps)     ! New Old-Infectious per day
        dimension nep(ips,jps)      ! New Exposed per day
        dimension ndh(ips,jps)      ! New Dead per day     
        dimension apr(ips,jps)      !Production of virus  
        dimension amvsc(ips,jps)    !Movement probability of Susceptible-Child  
        dimension amvsy(ips,jps)    !Movement probability of Susceptible-Young  
        dimension amvso(ips,jps)    !Movement probability of Susceptible-Old  
        dimension amvic(ips,jps)    !Movement probability of Infectious-Child  
        dimension amviy(ips,jps)    !Movement probability of Infectious-Young  
        dimension amvio(ips,jps)    !Movement probability of Infectious-Old  
        dimension amvec(ips,jps)    !Movement probability of Exposed-Child   
        dimension amvey(ips,jps)    !Movement probability of Exposed-Young 
        dimension amveo(ips,jps)    !Movement probability of Exposed-Old  
        dimension mnipa_data(5000)
        dimension mndha_data(5000)
        dimension drop(ips,jps)     !Death rate for old
        dimension dryp(ips,jps)     !Death rate for young
        dimension drc(ips,jps)      !Death rate for child
        dimension ath(ips,jps)      !Threshold level for Susceptible
        dimension atb(ips,jps)      !Threshold level for Recovered
        dimension itmrv(ips,jps)    !Time since Recovered
        dimension nv(ips,jps)       !Vaccinated 
        dimension ncv(ips,jps)      !Vaccinated-Child
        dimension nyv(ips,jps)      !Vaccinated-Young 
        dimension nov(ips,jps)      !Vaccinated-Old

c	idum=0.0
        call random_seed
        iseed=101  
        iwrt=0.0   
        max_ma=0
        max_ndh=0


        open(14,file='covid_UK_input.in',status='old')	
        open(151,file='new_case_p_d.dat',status='unknown')	! infection daily
        open(152,file='new_death_p_d.dat',status='unknown')	! dead daily
        open(251,file='total_cases.dat',status='unknown')	! cumulative infection
        open(99251,file='total_death.dat',status='unknown')	! cumulative dead
        open(151151,file='wdc.dat',status='unknown')		! child dead daily
        open(151152,file='wdy.dat',status='unknown')		! young dead daily
        open(151153,file='wdo.dat',status='unknown')		! old dead daily
        open(151154,file='infect_c.dat',status='unknown')	!child infected daily
        open(151155,file='infect_y.dat',status='unknown')	!young infected daily
        open(151156,file='infect_o.dat',status='unknown')	!old infected daily

c      Used parameter of our simulation
        read(14,*)dh		! Step length for spatial integration
        read(14,*)rad		! Radius for initial distribution of infectious
        read(14,*)dv		! Diffusion constant of virus
        read(14,*)phi		! Degradation rate of virus
        read(14,*)theta		! Production rate of virus
        read(14,*)kt,ktr	! total number of iteration
        read(14,*)dt        	! Time integration step length
        read(14,*)amv_sy        ! movement probability of young 
        read(14,*)amv_sc       	! movement probability of children
        read(14,*)amv_so       	! movement probability of old 
        read(14,*)amv_ec        ! movement probability of infected child
        read(14,*)amv_ey        ! movement probability of infected young
        read(14,*)amv_eo        ! movement probability of infected young 
        read(14,*)amv_ic        ! movement probability of infected child
        read(14,*)amv_iy        ! movement probability of infected young
        read(14,*)amv_io        ! movement probability of infected young
        read(14,*)flckdn1       ! Multiplication factor for lockdown-1
        read(14,*)flckdn2       ! Multiplication factor for lockdown-2
        read(14,*)flckdn3       ! Multiplication factor for lockdown-3
        read(14,*)hp        	! random number for healthy people(initialization)
        read(14,*)ainp        	! random number for infected people(initialization)
        read(14,*)thv        	! Threshold level of virus for infection
        read(14,*)it1d          ! Day-1
        read(14,*)it13d         ! Day-13
        read(14,*)it18d         ! Day-18
        read(14,*)it330d        ! Day-330
        read(14,*)it480d        ! Day-480
        read(14,*)it25d         ! Day-25
        read(14,*)it100d        ! Day-100
        read(14,*)it320d        ! Day-320
        read(14,*)it90d         ! Day-90
        read(14,*)it340d        ! Day-340
        read(14,*)it370d        ! Day-370
        read(14,*)it400d        ! Day-400
        read(14,*)it535d        ! Day-535
        read(14,*)it565d        ! Day-565
        read(14,*)it28d         ! Day-28
        read(14,*)it96d         ! Day-96
        read(14,*)it188d        ! Day-188
        read(14,*)it252d        ! Day-252
        read(14,*)it289d        ! Day-289
        read(14,*)it313d        ! Day-313
        read(14,*)it470d        ! Day-470
        read(14,*)precov        ! probability of recovery
        read(14,*)ainfc        	! infection rate of child
        read(14,*)ainfy        	! infection rate of young
        read(14,*)ainfo        	! infection rate of old
        read(14,*)it_incu       ! incubation period
        read(14,*)avca        	! probability of vaccination
        read(14,*)pcrt		    ! probability of critical illness

c    Initialization
      do i=1,ips
       do j=1,jps
        ns(i,j)=0.0
        nsc(i,j)=0.0
        nsy(i,j)=0.0
        nso(i,j)=0.0
        ne(i,j)=0.0
        ni(i,j)=0.0
        nic(i,j)=0.0
        niy(i,j)=0.0
        nio(i,j)=0.0
        nr(i,j)=0.0
        nrc(i,j)=0.0
        nry(i,j)=0.0
        nro(i,j)=0.0
        neo(i,j)=0.0
        ney(i,j)=0.0
        nec(i,j)=0.0
        vir0(i,j)=0.0
        itm(i,j)=0.0
        itmc(i,j)=0.0
        itmy(i,j)=0.0
        itmo(i,j)=0.0
        itmrv(i,j)=0.0
        dth(i,j)=0.0
        dthc(i,j)=0.0
        dthy(i,j)=0.0
        dtho(i,j)=0.0
        nip(i,j)=0.0
        nipc(i,j)=0.0
        nipy(i,j)=0.0
        nipo(i,j)=0.0
        nep(i,j)=0
        ndh(i,j)=0.0
        nrr(i,j)=0.0
        apr(i,j)=theta
        amvsc(i,j)=amv_sc
        amvsy(i,j)=amv_sy
        amvso(i,j)=amv_so
        amvic(i,j)=amv_ic
        amviy(i,j)=amv_iy
        amvio(i,j)=amv_io
        amvec(i,j)=amv_ec
        amvey(i,j)=amv_ey
        amveo(i,j)=amv_eo
        drop(i,j)=0.0
        dryp(i,j)=0.0
        drc(i,j)=0.0
        ath(i,j)=thv
        nov(i,j)=0.0
        nyv(i,j)=0.0
        nv(i,j)=0.0
       enddo
      enddo

      sumns=0
      sumne=0
      sumnc=0
      sumny=0
      sumno=0
      sumnec=0
      sumney=0
      sumneo=0
      ici=1

      do 23 i=1,ips
       do 24 j=1,jps
        x=i*dh
        y=j*dh
        r=sqrt((x-1.0)**2+(y-1.0)**2)
c      **************************************************************
       !Initialization of susceptible agents
c      **************************************************************
c      ar=rand(idum)
        call random_number(harvest=ar)
        if(ar.le.hp)then
         ns(i,j)=1.0
         sumns=sumns+1
        else
         ns(i,j)=0
        endif
c      **************************************************************
       !Initialization of Exposed agents
c      **************************************************************           
        if(ici.eq.1)then
         if(r.le.rad)then
c       arr=rand(idum)
          call random_number(harvest=arr)
          if(ns(i,j).eq.0.and.arr.le.ainp)then
           if(sumne.lt.3)then
            ne(i,j)=1.0
            sumne=sumne+1
           endif
          else
           ne(i,j)=0.0
          endif
         endif !
        endif

        if(sumne.eq.3)then
         ici=0
        endif
c      **************************************************************
        !Initialization of virus
        vir0(i,j)=0.0
c       **************************************************************
c        age distribution of Susceptible agents:25%-child,50%-young,25%-old
c       **************************************************************
        if(ns(i,j).eq.1)then
c       ar1=rand(idum)
         call random_number(harvest=ar1)
         if(ar1.le.0.25)then
          nsc(i,j)=1  !Child
          sumnc=sumnc+1
         else if(ar1.gt.0.25.and.ar1.lt.0.75)then
          nsy(i,j)=1  !young
          sumny=sumny+1
         else
          nso(i,j)=1  !old 
          sumno=sumno+1
         endif
        endif
c      **************************************************************
c        age distribution of Exposed agents:25%-child,25%-young,50%-old
c      **************************************************************
        if(ne(i,j).eq.1)then
c      aar1=rand(idum)
         call random_number(harvest=aar1)
         if(aar1.le.0.5)then
          neo(i,j)=1  !Exposed old 
          sumneo=sumneo+1
         else if(aar1.gt.0.5.and.aar1.lt.0.75)then
          ney(i,j)=1  !Exposed young
          sumney=sumney+1
         else
          nec(i,j)=1  !Exposed Child
          sumnec=sumnec+1
         endif
        endif
c     *****************************************************************
24     continue
23    continue
      
      write(*,*)'ns=',sumns
      write(*,*)'ne=',sumne
      write(*,*)'nc=',sumnc
      write(*,*)'ny=',sumny
      write(*,*)'no=',sumno
      write(*,*)'nec=',sumnec
      write(*,*)'ney=',sumney
      write(*,*)'neo=',sumneo
      write(*,*)'total time(hour)=',(kt*dt)
c     *****************************************************************
c       #############################################################
      nipa=0
      nipac=0
      nipay=0
      nipao=0
      ndha=0
      icount=0
      icot=0
      icount1=0
      iflag=1
      iflag1=1
      iflag2=1
      ntca=0
      ntdh=0
      ndthc=0
      ndthy=0
      ndtho=0
      ! Diffusion parameters
      vv=dt*dv/dh**2
        write(*,*)'tolerance=dv',1.0/vv
c       ***************************************************************
c	Time loop starts here
c       ***************************************************************
      do 15 k=1,kt   

c      Infection rates are higher for new variants         
       if(k.ge.it480d)then !Delta variant of virus
        pinfc=ainfc*1.4 
        pinfy=ainfy*1.4
        pinfo=ainfo*1.4
       else
        pinfc=ainfc
        pinfy=ainfy
        pinfo=ainfo
       endif 

c      defining vaccination of old and young	
       if(k.le.it340d)then
        avcao=0.0*avca
        avcay=0.0*avca
       elseif(k.gt.it340d.and.k.le.it370d)then
        avcao=1.0*avca
        avcay=0.5*avca 
       elseif(k.gt.it370d.and.k.le.it400d)then
        avcao=1.5*avca
        avcay=0.75*avca
       else
        avcao=1.5*avca
        avcay=1.0*avca
       endif

c      defining vaccination of child
       if(k.le.it535d)then
        avcac=0.0*avca
       elseif(k.gt.it535d.and.k.le.it565d)then
        avcac=0.5*avca
       else
	avcac=1.0*avca
       endif
         
c       ***************************************************************
c	Lattice loop starts here
c       ***************************************************************
       do 18 i=1,ips
        do 19 j=1,jps
c       **************************************************************
        ! Implementation of zero flux boundary condition
         if(i.eq.1.and.j.eq.1)then
          il=i
          ir=i+1
          jb=j
          jt=j+1
         else if(i.eq.1.and.j.gt.1.and.j.lt.jps)then
          il=i
          ir=i+1
          jb=j-1
          jt=j+1
         else if(i.eq.1.and.j.eq.jps)then
          il=i
          ir=i+1
          jb=j-1
          jt=j
         else if(i.gt.1.and.i.lt.ips.and.j.eq.1)then
          il=i-1
          ir=i+1
          jb=j
          jt=j+1
         else if(i.gt.1.and.i.lt.ips.and.j.eq.jps)then
          il=i-1
          ir=i+1
          jb=j-1
          jt=j
         else if(i.eq.ips.and.j.eq.1)then
          il=i-1
          ir=i
          jb=j
          jt=j+1
         else if(i.eq.ips.and.j.gt.1.and.j.lt.jps)then
          il=i-1
          ir=i
          jb=j-1
          jt=j+1
         else if(i.eq.ips.and.j.eq.jps)then
          il=i-1
          ir=i
          jb=j-1
          jt=j
         else
          il=i-1
          ir=i+1
          jb=j-1
          jt=j+1 
         endif
c       **************************************************************

c       **************************************************************
c       Production rate of virus for infected and exposed agents
         if(ne(i,j).eq.1.and.ni(i,j).eq.1)then
          apr(i,j)=theta !production rate of virus
         else if(ne(i,j).eq.1.and.ni(i,j).eq.0)then
          apr(i,j)=0.68*theta
         endif
c       **************************************************************

	 ! Equation for Virus 	
         vir(i,j)=vir0(i,j)*(1.0-4.0*vv-phi*dt)+vv*(vir0(ir,j)
     $     +vir0(il,j)+vir0(i,jt)+vir0(i,jb))+ne(i,j)*apr(i,j)*dt    

         vir0(i,j)=vir(i,j)  
         x=i*dh
         y=j*dh

c       ***************************************************************
c      Implementation of lock-down
         if(k.lt.it28d)then !28days (1-28)day Unrestricted
          amvsc(i,j)=amv_sc
          amvsy(i,j)=amv_sy
          amvso(i,j)=amv_so
          amvec(i,j)=amv_ec
          amvey(i,j)=amv_ey
          amveo(i,j)=amv_eo  
         else if(k.ge.it28d.and.k.lt.it96d)then !68days (29-96)day Lock-down-1
          amvsc(i,j)=flckdn1*amv_sc
          amvsy(i,j)=flckdn1*amv_sy
          amvso(i,j)=flckdn1*amv_so
          amvec(i,j)=flckdn1*amv_ec
          amvey(i,j)=flckdn1*amv_ey
          amveo(i,j)=flckdn1*amv_eo 
         else if(k.ge.it96d.and.k.lt.it188d)then !92 days (97-188)day Lock-down-2
          amvsc(i,j)=flckdn2*amv_sc
          amvsy(i,j)=flckdn2*amv_sy
          amvso(i,j)=flckdn2*amv_so
          amvec(i,j)=flckdn2*amv_ec
          amvey(i,j)=flckdn2*amv_ey
          amveo(i,j)=flckdn2*amv_eo  
         else if(k.ge.it188d.and.k.lt.it252d)then !37days (189-252)day Lock-down-3
          amvsc(i,j)=flckdn3*amv_sc
          amvsy(i,j)=flckdn3*amv_sy
          amvso(i,j)=flckdn3*amv_so
          amvec(i,j)=flckdn3*amv_ec
          amvey(i,j)=flckdn3*amv_ey
          amveo(i,j)=flckdn3*amv_eo
         else if(k.ge.it252d.and.k.lt.it289d)then !37days (253-289)day Lock-down-2
          amvsc(i,j)=flckdn2*amv_sc
          amvsy(i,j)=flckdn2*amv_sy
          amvso(i,j)=flckdn2*amv_so
          amvec(i,j)=flckdn2*amv_ec
          amvey(i,j)=flckdn2*amv_ey
          amveo(i,j)=flckdn2*amv_eo 
         else if(k.ge.it289d.and.k.lt.it313d)then !24days (290-313)day Unrestricted
          amvsc(i,j)=amv_sc
          amvsy(i,j)=amv_sy
          amvso(i,j)=amv_so
          amvec(i,j)=amv_ec
          amvey(i,j)=amv_ey
          amveo(i,j)=amv_eo
         else if(k.ge.it313d.and.k.lt.it470d)then !157days (314-470)day Lock-down-2
          amvsc(i,j)=flckdn2*amv_sc
          amvsy(i,j)=flckdn2*amv_sy
          amvso(i,j)=flckdn2*amv_so
          amvec(i,j)=flckdn2*amv_ec
          amvey(i,j)=flckdn2*amv_ey
          amveo(i,j)=flckdn2*amv_eo 
         else ! 155days(471-625)days Lock-down-3
          amvsc(i,j)=flckdn3*amv_sc
          amvsy(i,j)=flckdn3*amv_sy
          amvso(i,j)=flckdn3*amv_so
          amvec(i,j)=flckdn3*amv_ec
          amvey(i,j)=flckdn3*amv_ey
          amveo(i,j)=flckdn3*amv_eo
         endif !if(k.lt.1200)then
c       ***************************************************************
    
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     Spreading of Covid-19
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     ****************************************************************
c     SUSCEPTIBLE to EXPOSED 
c     ****************************************************************
       if(vir(i,j).ge.thv)then	      ! corrsing virus threshold
        if(ns(i,j).eq.1.and.nv(i,j).eq.0)then  ! for susceptible & unvaccinated agents
         if(nrr(i,j).eq.0.and.nr(i,j).eq.0)then ! who are not recovered from exposed or infection     
          if(nso(i,j).eq.1)then 	!exposed: old
           ne(i,j)=1
           neo(i,j)=1
           nso(i,j)=0
           ns(i,j)=0
          else if(nsy(i,j).eq.1)then 	!exposed: young
           ne(i,j)=1
           ney(i,j)=1
           nsy(i,j)=0
           ns(i,j)=0
          else if(nsc(i,j).eq.1)then 	!exposed: children
           ne(i,j)=1
           nec(i,j)=1
           nsc(i,j)=0
           ns(i,j)=0
          endif
         else if(nrr(i,j).eq.1)then 	!for agents recovered from exposed
          if(itmrv(i,j).gt.it90d)then
           if(nso(i,j).eq.1)then 	!exposed: old
            ne(i,j)=1
            neo(i,j)=1
            nso(i,j)=0
            ns(i,j)=0
            nep(i,j)=1
            nrr(i,j)=0
            itmrv(i,j)=0
           else if(nsy(i,j).eq.1)then 	!exposed: young
            ne(i,j)=1
            ney(i,j)=1
            nsy(i,j)=0
            ns(i,j)=0
            nep(i,j)=1
            nrr(i,j)=0
            itmrv(i,j)=0
           else if(nsc(i,j).eq.1)then 	!exposed: child
            ne(i,j)=1
            nec(i,j)=1
            nsc(i,j)=0
            ns(i,j)=0
            nep(i,j)=1
            nrr(i,j)=0
            itmrv(i,j)=0
           endif
          endif
         else if(nr(i,j).eq.1)then	!for agents recovered from infected
          if(itmrv(i,j).gt.it90d)then
           if(nso(i,j).eq.1)then 	!exposed: old
            ne(i,j)=1
            neo(i,j)=1
            nso(i,j)=0
            ns(i,j)=0
            nep(i,j)=1
            nr(i,j)=0
            itmrv(i,j)=0
           else if(nsy(i,j).eq.1)then 	!exposed: young
            ne(i,j)=1
            ney(i,j)=1
            nsy(i,j)=0
            ns(i,j)=0
            nep(i,j)=1
            nr(i,j)=0
            itmrv(i,j)=0
           else if(nsc(i,j).eq.1)then 	!exposed: child
            ne(i,j)=1
            nec(i,j)=1
            nsc(i,j)=0
            ns(i,j)=0
            nep(i,j)=1
            nr(i,j)=0
            itmrv(i,j)=0
           endif
          endif
         endif ! recovered nr and nrr
        endif ! healthy & unvaccinated
      endif   !thv=threshold level of virus
c     ****************************************************************
c	Update of time since exposed
c     ****************************************************************
      if(ne(i,j).eq.1)then 
       itm(i,j)=itm(i,j)+1  ! day counting from exposed
      endif
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     EXPOSED to INFECTED 
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if(ne(i,j).eq.1.and.ni(i,j).eq.0.and.itm(i,j).ge.it_incu)then !exposed & incubation period	
c       cnf=rand(idum)
       call random_number(harvest=cnf)
       if(nec(i,j).eq.1)then		!infected: child
        if(cnf.lt.pinfc)then
         ni(i,j)=1
         nic(i,j)=1
         nip(i,j)=1
         nipc(i,j)=1
        else
         ne(i,j)=0
         nec(i,j)=0
         nsc(i,j)=1
         ns(i,j)=1
         nrr(i,j)=1
        endif
       else if(ney(i,j).eq.1)then	!infected: young
        if(cnf.lt.pinfy)then
         ni(i,j)=1
         niy(i,j)=1
         nip(i,j)=1
         nipy(i,j)=1
        else
         ne(i,j)=0
         ney(i,j)=0
         nsy(i,j)=1
         ns(i,j)=1
         nrr(i,j)=1
        endif
       else if(neo(i,j).eq.1)then	!infected: old
        if(cnf.lt.pinfo)then
         ni(i,j)=1
         nio(i,j)=1
         nip(i,j)=1
         nipo(i,j)=1
        else
         ne(i,j)=0
         neo(i,j)=0
         nso(i,j)=1
         ns(i,j)=1
         nrr(i,j)=1
        endif
       endif 
      endif !exposed & incubation period
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     CALCULATION OF RANDOM MOVEMENT
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if(mod(k,1).eq.0)then !chances of movement in every iteration=1hour
c     ****************************************************************
c     Movement for susceptible agents
c     ****************************************************************
      if(ns(i,j).eq.1)then 	!checking the presence of susceptible agent
c      pt=rand(idum) !Random number for movement of(Right,Left,Top and bottom)
       call random_number(harvest=pt)
      if(pt.le.0.25)then !Right
       if(ns(ir,j).eq.0.and.ne(ir,j).eq.0)then
c        ar2=rand(idum)
         call random_number(harvest=ar2)
        if(nsy(i,j).eq.1.and.ar2.lt.amvsy(i,j))then
         nsy(i,j)=0
         ns(i,j)=0
         nsy(ir,j)=1
         ns(ir,j)=1
         nr(ir,j)=nr(i,j)
         nr(i,j)=0
         nry(ir,j)=nry(i,j)
         nry(i,j)=0
         nrr(ir,j)=nrr(i,j)
         nrr(i,j)=0
         itmrv(ir,j)=itmrv(i,j)
         itmrv(i,j)=0
         nyv(ir,j)=nyv(i,j)
         nyv(i,j)=0
         nv(ir,j)=nv(i,j)
         nv(i,j)=0         
        else if(nsc(i,j).eq.1.and.ar2.lt.amvsc(i,j))then
         nsc(i,j)=0
         ns(i,j)=0
         nsc(ir,j)=1
         ns(ir,j)=1
         nr(ir,j)=nr(i,j)
         nr(i,j)=0
         nrc(ir,j)=nrc(i,j)
         nrc(i,j)=0
         nrr(ir,j)=nrr(i,j)
         nrr(i,j)=0
         itmrv(ir,j)=itmrv(i,j)
         itmrv(i,j)=0
         ncv(ir,j)=ncv(i,j)
         ncv(i,j)=0
         nv(ir,j)=nv(i,j)
         nv(i,j)=0
        else if(nso(i,j).eq.1.and.ar2.lt.amvso(i,j))then
         nso(i,j)=0
         ns(i,j)=0
         nso(ir,j)=1
         ns(ir,j)=1
         nr(ir,j)=nr(i,j)
         nr(i,j)=0
         nro(ir,j)=nro(i,j)
         nro(i,j)=0
         nrr(ir,j)=nrr(i,j)
         nrr(i,j)=0
         itmrv(ir,j)=itmrv(i,j)
         itmrv(i,j)=0
         nov(ir,j)=nov(i,j)
         nov(i,j)=0
         nv(ir,j)=nv(i,j)
         nv(i,j)=0
        endif !if(nsy(i,j).eq.1.and.ar2.lt.amv_sy)then
       endif !if(ns(ir,j).eq.0.and.ne(ir,j).eq.0)then

      else if(pt.gt.0.25.and.pt.le.0.5)then !Left
       if(ns(il,j).eq.0.and.ne(il,j).eq.0)then
c        ar3=rand(idum)
         call random_number(harvest=ar3)
        if(nsy(i,j).eq.1.and.ar3.lt.amvsy(i,j))then
         nsy(i,j)=0
         ns(i,j)=0
         nsy(il,j)=1
         ns(il,j)=1
         nr(il,j)=nr(i,j)
         nr(i,j)=0
         nry(il,j)=nry(i,j)
         nry(i,j)=0
         nrr(il,j)=nrr(i,j)
         nrr(i,j)=0
         itmrv(il,j)=itmrv(i,j)
         itmrv(i,j)=0
         nyv(il,j)=nyv(i,j)
         nyv(i,j)=0
         nv(il,j)=nv(i,j)
         nv(i,j)=0         
        else if(nsc(i,j).eq.1.and.ar3.lt.amvsc(i,j))then
         nsc(i,j)=0
         ns(i,j)=0
         nsc(il,j)=1
         ns(il,j)=1
         nr(il,j)=nr(i,j)
         nr(i,j)=0
         nrc(il,j)=nrc(i,j)
         nrc(i,j)=0
         nrr(il,j)=nrr(i,j)
         nrr(i,j)=0
         itmrv(il,j)=itmrv(i,j)
         itmrv(i,j)=0
         ncv(il,j)=ncv(i,j)
         ncv(i,j)=0
         nv(il,j)=nv(i,j)
         nv(i,j)=0
        else if(nso(i,j).eq.1.and.ar3.lt.amvso(i,j))then
         nso(i,j)=0
         ns(i,j)=0
         nso(il,j)=1
         ns(il,j)=1
         nr(il,j)=nr(i,j)
         nr(i,j)=0
         nro(il,j)=nro(i,j)
         nro(i,j)=0
         nrr(il,j)=nrr(i,j)
         nrr(i,j)=0
         itmrv(il,j)=itmrv(i,j)
         itmrv(i,j)=0
         nov(il,j)=nov(i,j)
         nov(i,j)=0
         nv(il,j)=nv(i,j)
         nv(i,j)=0        
        endif !if(nsy(i,j).eq.1.and.ar3.lt.amv_sy)then
       endif !if(ns(il,j).eq.0.and.ne(il,j).eq.0)then

      else if(pt.gt.0.5.and.pt.le.0.75)then !Top
       if(ns(i,jt).eq.0.and.ne(i,jt).eq.0)then
c        ar4=rand(idum)
              call random_number(harvest=ar4)
        if(nsy(i,j).eq.1.and.ar4.lt.amvsy(i,j))then
         nsy(i,j)=0
         ns(i,j)=0
         nsy(i,jt)=1
         ns(i,jt)=1
         nr(i,jt)=nr(i,j)
         nr(i,j)=0
         nry(i,jt)=nry(i,j)
         nry(i,j)=0
         nrr(i,jt)=nrr(i,j)
         nrr(i,j)=0
         itmrv(i,jt)=itmrv(i,j)
         itmrv(i,j)=0
         nyv(i,jt)=nyv(i,j)
         nyv(i,j)=0
         nv(i,jt)=nv(i,j)
         nv(i,j)=0         
        else if(nsc(i,j).eq.1.and.ar4.lt.amvsc(i,j))then
         nsc(i,j)=0
         ns(i,j)=0
         nsc(i,jt)=1
         ns(i,jt)=1
         nr(i,jt)=nr(i,j)
         nr(i,j)=0
         nrc(i,jt)=nrc(i,j)
         nrc(i,j)=0
         nrr(i,jt)=nrr(i,j)
         nrr(i,j)=0
         itmrv(i,jt)=itmrv(i,j)
         itmrv(i,j)=0
         ncv(i,jt)=ncv(i,j)
         ncv(i,j)=0
         nv(i,jt)=nv(i,j)
         nv(i,j)=0
        else if(nso(i,j).eq.1.and.ar4.lt.amvso(i,j))then
         nso(i,j)=0
         ns(i,j)=0
         nso(i,jt)=1
         ns(i,jt)=1
         nr(i,jt)=nr(i,j)
         nr(i,j)=0
         nro(i,jt)=nro(i,j)
         nro(i,j)=0
         nrr(i,jt)=nrr(i,j)
         nrr(i,j)=0
         itmrv(i,jt)=itmrv(i,j)
         itmrv(i,j)=0
         nov(i,jt)=nov(i,j)
         nov(i,j)=0
         nv(i,jt)=nv(i,j)
         nv(i,j)=0         
        endif !if(nsy(i,j).eq.1.and.ar4.lt.amv_sy)then
       endif !if(ns(i,jt).eq.0.and.ne(i,jt).eq.0)then

      else if(pt.gt.0.75)then !Bottom
       if(ns(i,jb).eq.0.and.ne(i,jb).eq.0)then
c        ar5=rand(idum)
         call random_number(harvest=ar5)
        if(nsy(i,j).eq.1.and.ar5.lt.amvsy(i,j))then
         nsy(i,j)=0
         ns(i,j)=0
         nsy(i,jb)=1
         ns(i,jb)=1
         nr(i,jb)=nr(i,j)
         nr(i,j)=0
         nry(i,jb)=nry(i,j)
         nry(i,j)=0
         nrr(i,jb)=nrr(i,j)
         nrr(i,j)=0
         itmrv(i,jb)=itmrv(i,j)
         itmrv(i,j)=0
         nyv(i,jb)=nyv(i,j)
         nyv(i,j)=0
         nv(i,jb)=nv(i,j)
         nv(i,j)=0         
        else if(nsc(i,j).eq.1.and.ar5.lt.amvsc(i,j))then
         nsc(i,j)=0
         ns(i,j)=0
         nsc(i,jb)=1
         ns(i,jb)=1
         nr(i,jb)=nr(i,j)
         nr(i,j)=0
         nrc(i,jb)=nrc(i,j)
         nrc(i,j)=0
         nrr(i,jb)=nrr(i,j)
         nrr(i,j)=0
         itmrv(i,jb)=itmrv(i,j)
         itmrv(i,j)=0
         ncv(i,jb)=ncv(i,j)
         ncv(i,j)=0
         nv(i,jb)=nv(i,j)
         nv(i,j)=0
        else if(nso(i,j).eq.1.and.ar5.lt.amvso(i,j))then
         nso(i,j)=0
         ns(i,j)=0
         nso(i,jb)=1
         ns(i,jb)=1
         nr(i,jb)=nr(i,j)
         nr(i,j)=0
         nro(i,jb)=nro(i,j)
         nro(i,j)=0
         nrr(i,jb)=nrr(i,j)
         nrr(i,j)=0
         itmrv(i,jb)=itmrv(i,j)
         itmrv(i,j)=0
         nov(i,jb)=nov(i,j)
         nov(i,j)=0
         nv(i,jb)=nv(i,j)
         nv(i,j)=0         
        endif !if(nsy(i,j).eq.1.and.ar4.lt.amv_sy)then
       endif !if(ns(i,jt).eq.0.and.ne(i,jt).eq.0)then

      endif !if(pt.le.0.25)then !Right  
c    *****************************************************************
      endif !if(ns(i,j).eq.1.and.nr(i,j).eq.0)then 
c    *****************************************************************
c    Movement for infected and exposed agents
c     ****************************************************************
      if(ne(i,j).eq.1)then !checking the presence of infected or exposed people
c    *****************************************************************
c      pt1=rand(idum) !Random number for movement of(Right,Left,Top and bottom)
       call random_number(harvest=pt1)
      if(pt1.le.0.25)then !Right
       if(ns(ir,j).eq.0.and.ne(ir,j).eq.0)then
c        arr2=rand(idum)
      call random_number(harvest=arr2)
      if(nec(i,j).eq.1.and.ni(i,j).eq.0.and.arr2.lt.amvec(i,j))then!exposed child
         ne(i,j)=0
         nec(i,j)=0
         ne(ir,j)=1
         nec(ir,j)=1
         itm(ir,j)=itm(i,j)
         itm(i,j)=0
      else if(nec(i,j).eq.1.and.ni(i,j).eq.1.and.arr2.lt.amvic(i,j))then!infected child
         ne(i,j)=0
         nec(i,j)=0
         ne(ir,j)=1
         nec(ir,j)=1
         itm(ir,j)=itm(i,j)
         itm(i,j)=0
         ni(ir,j)=1
         ni(i,j)=0
         nic(ir,j)=1
         nic(i,j)=0
      endif

      if(ney(i,j).eq.1.and.ni(i,j).eq.0.and.arr2.lt.amvey(i,j))then!exposed young
         ne(i,j)=0
         ney(i,j)=0
         ne(ir,j)=1
         ney(ir,j)=1
         itm(ir,j)=itm(i,j)
         itm(i,j)=0
      else if(ney(i,j).eq.1.and.ni(i,j).eq.1.and.arr2.lt.amviy(i,j))then!infected young
         ne(i,j)=0
         ney(i,j)=0
         ne(ir,j)=1
         ney(ir,j)=1
         itm(ir,j)=itm(i,j)
         itm(i,j)=0
         ni(ir,j)=1
         ni(i,j)=0
         niy(ir,j)=1
         niy(i,j)=0
      endif

      if(neo(i,j).eq.1.and.ni(i,j).eq.0.and.arr2.lt.amveo(i,j))then!exposed old
         ne(i,j)=0
         neo(i,j)=0
         ne(ir,j)=1
         neo(ir,j)=1
         itm(ir,j)=itm(i,j)
         itm(i,j)=0
      else if(neo(i,j).eq.1.and.ni(i,j).eq.1.and.arr2.lt.amvio(i,j))then!infected old
         ne(i,j)=0
         neo(i,j)=0
         ne(ir,j)=1
         neo(ir,j)=1
         itm(ir,j)=itm(i,j)
         itm(i,j)=0
         ni(ir,j)=1
         ni(i,j)=0
         nio(ir,j)=1
         nio(i,j)=0
      endif
       endif !if(ns(ir,j).eq.0.and.ne(ir,j).eq.0)then

      else if(pt1.gt.0.25.and.pt1.le.0.5)then !Left
       if(ns(il,j).eq.0.and.ne(il,j).eq.0)then
c        arr3=rand(idum)
         call random_number(harvest=arr3)
      if(nec(i,j).eq.1.and.ni(i,j).eq.0.and.arr3.lt.amvec(i,j))then
         ne(i,j)=0
         nec(i,j)=0
         ne(il,j)=1
         nec(il,j)=1
         itm(il,j)=itm(i,j)
         itm(i,j)=0
      else if(nec(i,j).eq.1.and.ni(i,j).eq.1.and.arr3.lt.amvic(i,j))then
         ne(i,j)=0
         nec(i,j)=0
         ne(il,j)=1
         nec(il,j)=1
         itm(il,j)=itm(i,j)
         itm(i,j)=0 
         ni(il,j)=1
         ni(i,j)=0
         nic(il,j)=1
         nic(i,j)=0
      endif

      if(ney(i,j).eq.1.and.ni(i,j).eq.0.and.arr3.lt.amvey(i,j))then
         ne(i,j)=0
         ney(i,j)=0
         ne(il,j)=1
         ney(il,j)=1
         itm(il,j)=itm(i,j)
         itm(i,j)=0
      else if(ney(i,j).eq.1.and.ni(i,j).eq.1.and.arr3.lt.amviy(i,j))then
         ne(i,j)=0
         ney(i,j)=0
         ne(il,j)=1
         ney(il,j)=1
         itm(il,j)=itm(i,j)
         itm(i,j)=0
         ni(il,j)=1
         ni(i,j)=0
         niy(il,j)=1
         niy(i,j)=0
      endif

      if(neo(i,j).eq.1.and.ni(i,j).eq.0.and.arr3.lt.amveo(i,j))then
         ne(i,j)=0
         neo(i,j)=0
         ne(il,j)=1
         neo(il,j)=1
         itm(il,j)=itm(i,j)
         itm(i,j)=0
      else if(neo(i,j).eq.1.and.ni(i,j).eq.1.and.arr3.lt.amvio(i,j))then
         ne(i,j)=0
         neo(i,j)=0
         ne(il,j)=1
         neo(il,j)=1
         itm(il,j)=itm(i,j)
         itm(i,j)=0
         ni(il,j)=1
         ni(i,j)=0
         nio(il,j)=1
         nio(i,j)=0
        endif 
       endif !if(ns(il,j).eq.0.and.ne(il,j).eq.0)then

      else if(pt1.gt.0.5.and.pt1.le.0.75)then !Top
       if(ns(i,jt).eq.0.and.ne(i,jt).eq.0)then
c        arr4=rand(idum)
      call random_number(harvest=arr4)
      if(nec(i,j).eq.1.and.ni(i,j).eq.0.and.arr4.lt.amvec(i,j))then
         ne(i,j)=0
         nec(i,j)=0
         ne(i,jt)=1
         nec(i,jt)=1
         itm(i,jt)=itm(i,j)
         itm(i,j)=0
      else if(nec(i,j).eq.1.and.ni(i,j).eq.1.and.arr4.lt.amvic(i,j))then
         ne(i,j)=0
         nec(i,j)=0
         ne(i,jt)=1
         nec(i,jt)=1
         itm(i,jt)=itm(i,j)
         itm(i,j)=0
         ni(i,jt)=1
         ni(i,j)=0
         nic(i,jt)=1
         nic(i,j)=0
      endif 

      if(ney(i,j).eq.1.and.ni(i,j).eq.0.and.arr4.lt.amvey(i,j))then
         ne(i,j)=0
         ney(i,j)=0
         ne(i,jt)=1
         ney(i,jt)=1
         itm(i,jt)=itm(i,j)
         itm(i,j)=0
      else if(ney(i,j).eq.1.and.ni(i,j).eq.1.and.arr4.lt.amviy(i,j))then
         ne(i,j)=0
         ney(i,j)=0
         ne(i,jt)=1
         ney(i,jt)=1
         itm(i,jt)=itm(i,j)
         itm(i,j)=0
         ni(i,jt)=1
         ni(i,j)=0
         niy(i,jt)=1
         niy(i,j)=0
      endif

      if(neo(i,j).eq.1.and.ni(i,j).eq.0.and.arr4.lt.amveo(i,j))then
         ne(i,j)=0
         neo(i,j)=0
         ne(i,jt)=1
         neo(i,jt)=1
         itm(i,jt)=itm(i,j)
         itm(i,j)=0
      else if(neo(i,j).eq.1.and.ni(i,j).eq.1.and.arr4.lt.amvio(i,j))then
         ne(i,j)=0
         neo(i,j)=0
         ne(i,jt)=1
         neo(i,jt)=1
         itm(i,jt)=itm(i,j)
         itm(i,j)=0
         ni(i,jt)=1
         ni(i,j)=0
         nio(i,jt)=1
         nio(i,j)=0
        endif
       endif !if(ns(i,jt).eq.0.and.ne(i,jt).eq.0)then

      else if(pt1.gt.0.75)then !Bottom
       if(ns(i,jb).eq.0.and.ne(i,jb).eq.0)then
c        arr5=rand(idum)
      call random_number(harvest=arr5)
      if(nec(i,j).eq.1.and.ni(i,j).eq.0.and.arr5.lt.amvec(i,j))then
         ne(i,j)=0
         nec(i,j)=0
         ne(i,jb)=1
         nec(i,jb)=1
         itm(i,jb)=itm(i,j)
         itm(i,j)=0
      else if(nec(i,j).eq.1.and.ni(i,j).eq.1.and.arr5.lt.amvic(i,j))then
         ne(i,j)=0
         nec(i,j)=0
         ne(i,jb)=1
         nec(i,jb)=1
         itm(i,jb)=itm(i,j)
         itm(i,j)=0
         ni(i,jb)=1
         ni(i,j)=0
         nic(i,jb)=1
         nic(i,j)=0
      endif

      if(ney(i,j).eq.1.and.ni(i,j).eq.0.and.arr5.lt.amvey(i,j))then
         ne(i,j)=0
         ney(i,j)=0
         ne(i,jb)=1
         ney(i,jb)=1
         itm(i,jb)=itm(i,j)
         itm(i,j)=0
      else if(ney(i,j).eq.1.and.ni(i,j).eq.1.and.arr5.lt.amviy(i,j))then
         ne(i,j)=0
         ney(i,j)=0
         ne(i,jb)=1
         ney(i,jb)=1
         itm(i,jb)=itm(i,j)
         itm(i,j)=0
         ni(i,jb)=1
         ni(i,j)=0
         niy(i,jb)=1
         niy(i,j)=0
      endif

      if(neo(i,j).eq.1.and.ni(i,j).eq.0.and.arr5.lt.amveo(i,j))then
         ne(i,j)=0
         neo(i,j)=0
         ne(i,jb)=1
         neo(i,jb)=1
         itm(i,jb)=itm(i,j)
         itm(i,j)=0
      else if(neo(i,j).eq.1.and.ni(i,j).eq.1.and.arr5.lt.amvio(i,j))then
         ne(i,j)=0
         neo(i,j)=0
         ne(i,jb)=1
         neo(i,jb)=1
         itm(i,jb)=itm(i,j)
         itm(i,j)=0
         ni(i,jb)=1
         ni(i,j)=0
         nio(i,jb)=1
         nio(i,j)=0
        endif !if(nec(i,j).eq.1.and.arr5.lt.amv_ec)then
       endif !if(ns(i,jb).eq.0.and.ne(i,jb).eq.0)then

      endif !if(pt1.le.0.25)then !Right   
c    *****************************************************************
      endif !if(ne(i,j).eq.1)then !checking the presence of infected people
c    *********************************************************************
c    *********************************************************************
      endif !if(mod(k,1).eq.0)then !chances of movement in every iteration=1hour
c    *********************************************************************
c     MOVEMENT CALCUALTION ENDS HERE
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     RECOVERY FROM INFECTED
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if(itm(i,j).ge.it13d.and.itm(i,j).le.it18d)then !recovery btwn 13 and 18 day
c       arv1=rand(idum)
       call random_number(harvest=arv1)
       if(arv1.le.precov)then
        if(neo(i,j).eq.1)then !recovery rate of old people is 
         ns(i,j)=1
         nso(i,j)=1
         nr(i,j)=1
         nro(i,j)=1
         ne(i,j)=0
         neo(i,j)=0
         itm(i,j)=0
         ni(i,j)=0
         nio(i,j)=0
        else if(ney(i,j).eq.1)then !recovery rate of young people
         ns(i,j)=1
         nsy(i,j)=1
         nr(i,j)=1
         nry(i,j)=1
         ne(i,j)=0
         ney(i,j)=0
         itm(i,j)=0
         ni(i,j)=0
         niy(i,j)=0
        else if(nec(i,j).eq.1)then !recovery rate of children
         ns(i,j)=1
         nsc(i,j)=1
         nr(i,j)=1
         nrc(i,j)=1
         ne(i,j)=0
         nec(i,j)=0
         itm(i,j)=0
         ni(i,j)=0
         nic(i,j)=0
        endif
       else
        go to 2
       endif
      endif
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
2	    ext=1 
c    *****************************************************************
c     Post recovery age(time) for the infected and exposed recovered
c    ***************************************************************** 
      if(nr(i,j).eq.1.or.nrr(i,j).eq.1)then 
         itmrv(i,j)=itmrv(i,j)+1
      endif
c    *****************************************************************

c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     HOSPITALIZED to DEAD
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if(ne(i,j).eq.1.and.itm(i,j).gt.it25d)then
c       arv1=rand(idum)
       call random_number(harvest=arv1)
       if(arv1.le.pcrt)then !60% chance of critical illness
        if(k.le.it100d)then
         drop(i,j)=0.90
         dryp(i,j)=0.70
         drc(i,j)=0.30
        elseif(k.gt.it100d.and.k.le.it320d)then
         drop(i,j)=0.70
         dryp(i,j)=0.50
         drc(i,j)=0.10
        else
         drop(i,j)=0.20
         dryp(i,j)=0.10
         drc(i,j)=0.05
        endif
c        dthr=rand(idum)
        call random_number(harvest=dthr)
        if(neo(i,j).eq.1.and.dthr.le.drop(i,j))then 
           ne(i,j)=0
           neo(i,j)=0
           dth(i,j)=1
           dtho(i,j)=1
           itmo(i,j)=itm(i,j)
           itm(i,j)=0
           ni(i,j)=0
           nio(i,j)=0
           ndh(i,j)=1

          else if(ney(i,j).eq.1.and.dthr.le.dryp(i,j))then 
           ne(i,j)=0
           ney(i,j)=0
           dth(i,j)=1
           dthy(i,j)=1
           itmy(i,j)=itm(i,j)
           itm(i,j)=0
           ni(i,j)=0
           niy(i,j)=0
           ndh(i,j)=1

          else if(nec(i,j).eq.1.and.dthr.le.drc(i,j))then 
           ne(i,j)=0
           nec(i,j)=0
           dth(i,j)=1
           dthc(i,j)=1
           itmc(i,j)=itm(i,j)
           itm(i,j)=0
           ni(i,j)=0
           nic(i,j)=0
           ndh(i,j)=1
          endif

         else 		! (1-prct) chance of recovery from critical illness
          if(neo(i,j).eq.1)then !recovery rate of old people is 
           ns(i,j)=1
           nso(i,j)=1
           nr(i,j)=1
           nro(i,j)=1
           ne(i,j)=0
           neo(i,j)=0
           itm(i,j)=0
           ni(i,j)=0
           nio(i,j)=0
          else if(ney(i,j).eq.1)then !recovery rate of young people
           ns(i,j)=1
           nsy(i,j)=1
           nr(i,j)=1
           nry(i,j)=1
           ne(i,j)=0
           ney(i,j)=0
           itm(i,j)=0
           ni(i,j)=0
           niy(i,j)=0
          else if(nec(i,j).eq.1)then !recovery rate of children
           ns(i,j)=1
           nsc(i,j)=1
           nr(i,j)=1
           nrc(i,j)=1
           ne(i,j)=0
           nec(i,j)=0
           itm(i,j)=0
           ni(i,j)=0
           nic(i,j)=0
          endif
         endif  !if(arv1.le.0.5)then

       endif !if(ne(i,j).eq.1.and.itm(i,j).gt.it25d)then
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c    *********************************************************************
c     VACCINATION  
c    *********************************************************************
       if(k.ge.it340d.and.mod(k,it1d).eq.0)then
       if(ns(i,j).eq.1)then
        call random_number(harvest=avc)           
         if(nso(i,j).eq.1.and.nov(i,j).eq.0.and.avc.le.avcao)then
          nov(i,j)=1
          nv(i,j)=1
         elseif(nsy(i,j).eq.1.and.nyv(i,j).eq.0.and.avc.le.avcay)then
          nyv(i,j)=1
          nv(i,j)=1
         endif
       endif!if(ns(i,j).eq.1)then
       endif!if(k.ge.it340d.and.mod(k,it1d).eq.0)then

       if(ns(i,j).eq.1.and.nsc(i,j).eq.1)then
       if(k.ge.it535d.and.mod(k,it1d).eq.0)then
         if(avc.le.avca)then
          ncv(i,j)=1
          nv(i,j)=1
         endif
       endif
       endif

c    &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
19        continue
18      continue
c      ****************************************************************
c      Data collection with time
c      ****************************************************************
      tm=(k*dt/it1d)
      nvir=0.0
      nenf=0
      nenfc=0
      nenfy=0
      nenfo=0
      nhelthy=0
      nhelthyc=0
      nhelthyy=0
      nhelthyo=0
      ndth=0
      nra=0
      nrca=0
      nrya=0
      nroa=0
      nia=0
      nica=0
      niya=0
      nioa=0
      nepa=0
      do i=1,ips
       do j=1,jps
        nvir=nvir+vir(i,j)         !virus
        nenf=nenf+ne(i,j)          !Exposed agents
        nenfc=nenfc+nec(i,j)       !Exposed child
        nenfy=nenfy+ney(i,j)       !Exposed young 
        nenfo=nenfo+neo(i,j)       !Exposed old 
        nhelthy=nhelthy+ns(i,j)    !Susceptible agents
        nhelthyc=nhelthyc+nsc(i,j) !Susceptible child
        nhelthyy=nhelthyy+nsy(i,j) !Susceptible young 
        nhelthyo=nhelthyo+nso(i,j) !Susceptible old
        ndth=ndth+dth(i,j)         !Dead agents
        ndthc=ndthc+dthc(i,j)      !Dead child
        ndthy=ndthy+dthy(i,j)      !Dead young
        ndtho=ndtho+dtho(i,j)      !Dead old
        nra=nra+nr(i,j)            !Recovered agents 
        nrca=nrca+nrc(i,j)         !Recovered child
        nrya=nrya+nry(i,j)         !Recovered young
        nroa=nroa+nro(i,j)         !Recovered old
        nia=nia+ni(i,j)            !Infected agents
        nica=nica+nic(i,j)         !Infected child
        niya=niya+niy(i,j)         !Infected young
        nioa=nioa+nio(i,j)         !Infected old
        nepa=nepa+nep(i,j)         !Daily exposed
        nipa=nipa+nip(i,j)         !Daily infected
        nipac=nipac+nipc(i,j)      !Daily infected child
        nipay=nipay+nipy(i,j)      !Daily infected young
        nipao=nipao+nipo(i,j)      !Daily infected old
        ndha=ndha+ndh(i,j)         !Daily dead
       enddo
      enddo
c      ****************************************************************
c      Daily data collection
c      ****************************************************************
        if(mod(k,it1d).eq.0)then
          icount=icount+1
          nipa_Data=nipa
          ndha_Data=ndha
        write(151,*)icount,nipa_Data
        write(152,*)icount,ndha_Data
          ntca=ntca+nipa_Data
          ntdh=ntdh+ndha_Data
        write(251,*)icount,ntca
        write(99251,*)icount,ntdh
          mnipa_data(icount)=nipa_Data 
          mndha_data(icount)=ndha_Data 
          if(nipa_Data.gt.max_ma)then
            max_ma=nipa_Data
          endif
          if(ndha_Data.gt.max_ndh)then
            max_ndh=ndha_Data
          endif
         write(153,*)icount,max_ma  
         write(1531,*)icount,max_ndh  
          nipa=0
          ndha=0
        endif
       
        if(mod(k,it1d).eq.0)then
         icot=icot+1
         nch_data=ndthc
         nyn_data=ndthy
         nod_data=ndtho
         nic_data=nipac
         niy_data=nipay
         nio_data=nipao
         write(151151,*)icot,nch_data
         write(151152,*)icot,nyn_data
         write(151153,*)icot,nod_data
         write(151154,*)icot,nic_data
         write(151155,*)icot,niy_data
         write(151156,*)icot,nio_data
         ndthc=0
         ndthy=0
         ndtho=0
         nipac=0
         nipay=0
         nipao=0
        endif
  
c      ****************************************************************
c        Data collection with different time interval
c      ****************************************************************     
       if(mod(k,ktr).eq.0)then
         iwrt=iwrt+1
         do i=1,ips
          do j=1,jps
           nwr1=200+(iwrt-1)*100
           nwr2=201+(iwrt-1)*100
           nwr3=202+(iwrt-1)*100
           nwr4=203+(iwrt-1)*100
           nwr5=204+(iwrt-1)*100
           nwr6=205+(iwrt-1)*100
           write(nwr1,*)i*dh,j*dh,vir(i,j)
           write(nwr2,*)i*dh,j*dh,dth(i,j)
           write(nwr3,*)i*dh,j*dh,ni(i,j)
           write(nwr4,*)i*dh,j*dh,nr(i,j)
           write(nwr5,*)i*dh,j*dh,ns(i,j)
           write(nwr6,*)i*dh,j*dh,ne(i,j)-ni(i,j)
          enddo
          write(nwr1,*)
          write(nwr2,*)
          write(nwr3,*)
          write(nwr4,*)
          write(nwr5,*)
          write(nwr6,*)
         enddo
        endif
c      ****************************************************************
c      Re-initialization for perday data collection
c      ****************************************************************
        do i=1,ips
         do j=1,jps
          nep(i,j)=0.0
          nip(i,j)=0.0
          ndh(i,j)=0.0
          dtho(i,j)=0.0
          dthy(i,j)=0.0
          dthc(i,j)=0.0
          nipc(i,j)=0.0
          nipy(i,j)=0.0
          nipo(i,j)=0.0
         enddo
        enddo     
c       ***************************************************************
15      continue     ! time loop ends here
c      ****************************************************************
       stop
       end

