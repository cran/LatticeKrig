       subroutine dfind1d( x1, n1, x2, n2, delta2, ind, rd, Nmax, iflag)
       integer n1,n2,ind(Nmax,2)
       integer kk, i,j, ic
       real*8 x1(n1), x2(n2), delta2(n2), rd(Nmax), dtemp
c****   counter  for accumulating close points
        kk=0 
          do  15 i= 1, n1
            kksave= kk
            do 10 j =1,n2
c**** accumulate squared differences
              dtemp= (x1(i) - x2(j))**2 
                if( dtemp.gt.delta2(j)) goto 10
c****       dtemp is less than delta2 so save it as a close point
              kk=kk+1
c**** check if there is still array space 
              if( kk .gt. Nmax) then 
                iflag= -1
                return
              else
                ind(kk,1)= i
                ind(kk,2)= j
                rd(kk)= sqrt( dtemp)
              endif     
 10        continue          
 15      continue
         iflag=1
         Nmax=kk  
 20      continue
      end

  
