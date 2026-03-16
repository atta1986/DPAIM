module polar_interp
   implicit none
contains

   subroutine polar_weighted_avg(Vt3, Vt5, d, vvol3, vte3, vvol5, vte5, pol, polq, c11, c12, c22, aiw_out, aiwq_out)
      implicit none
      real*8,    intent(in)  :: Vt3,Vt5
      real,    intent(in)  :: d(:)
      real,    intent(in)  :: vvol3(:), vte3(:)
      real,    intent(in)  :: vvol5(:), vte5(:)
      real,    intent(in)  :: pol(:,:), polq(:,:)
      real*8,    intent(in)  :: c11, c12, c22
      real,    intent(out) :: aiw_out(23), aiwq_out(23)

      integer :: n, ia, ib, j
      real    :: da, db, w3, w5, w_3, w_5
      real    :: Va3, Vb3, Va5, Vb5
      real    :: s1a, s1b, r3a, r3b
      real    :: s2a, s2b, r5a, r5b
      real    :: aiw_a, aiw_b, aiwq_a, aiwq_b
      real    :: denom_3,denom_5, denom_3a, denom_3b, denom_5a, denom_5b, denom, delta

      !print*,'c11,c12,c22',c11, c12, c22
      n = size(d)
      
      delta = 1.0e-6
      
      ! --- pick closest (ia) and 2nd closest (ib) by d ---
      ia = 1
      do j = 2, n
         if (d(j) < d(ia)) ia = j
      end do

      ib = 1
      if (ib == ia) ib = 2
      do j = 1, n
         if (j /= ia) then
            if (d(j) < d(ib)) ib = j
         end if
      end do

      Va3 = vte3(ia); Vb3 = vte3(ib)
      Va5 = vte5(ia); Vb5 = vte5(ib)

      denom_3 = (Vb3 - Va3)
      denom_5 = (Vb5 - Va5)
      
      !denom_3a = abs(Vt3-Va3) + delta
      !denom_3b = abs(Vt3-Vb3) + delta
      !denom_5a = abs(Vt5-Va5) + delta
      !denom_5b = abs(Vt5-Vb5) + delta

      if (abs(denom_3)< 1.0e-12.or.abs(denom_5)< 1.0e-12) then
         w3 = 0.0
         w5 = 0.0
      elseif (Vt3 >= min(Va3,Vb3) .and. Vt3 <= max(Va3,Vb3)) then
         w3 = (Vt3 - Va3)/denom_3
         w5 = (Vt5 - Va5)/denom_5
     !    w3 = max(0.0, min(1.0, abs(w_3)))
     !    w5 = max(0.0, min(1.0, abs(w_5)))
     !    write(5,*)'w3,w5,denom_3,denom_5,Va3,Va5,Vb3,Vb5,Vt3,Vt5,vvol3(ia),vvol3(ib)',w3,w5,denom_3,denom_5,Va3,Va5,Vb3,Vb5,Vt3,& 
     !        Vt5,vvol3(ia),vvol3(ib)
     !    print*,'w3,w5'
      else
         w3 = 0.0
         w5 = 0.0
      end if

      ! --- compute scaled aiw/aiwq for both references then interpolate
      ! ---
      !write(5,*)'ia,ib,Va3,Vb3,Vt3',ia,ib,Va3,Vb3,Vt3
      !write(5,*)'w3,w5,denom_3,denom_5,Va3,Vb3,Vt3,Vt5,Va5,Vb5',w3,w5,denom_3,denom_5,Va3,Vb3,Vt3,Vt5,Va5,Vb5
      s1a = Vt3 / vte3(ia);  r3a = Vt3 / vvol3(ia)
      s1b = Vt3 / vte3(ib);  r3b = Vt3 / vvol3(ib)

      s2a = Vt5 / vte5(ia);  r5a = Vt5 / vvol5(ia)
      s2b = Vt5 / vte5(ib);  r5b = Vt5 / vvol5(ib)
      !write(5,*)'s1a,s1b,s2a,s2b',s1a,s1b,s2a,s2b
      !write(5,*)'r3a,r3b,r5a,r5b,pol(ia,1),pol(ib,1)',r3a,r3b,r5a,r5b,pol(ia,1),pol(ib,1)
      !write(5,*)'diff',(abs(1.0 - c11*s1b)*r3b*pol(ib,1))-(abs(1.0 - c11*s1a)*r3a*pol(ia,1))
      !write(5,*)'ia,ib,s2a,r5a,polq(2,1),aiwq_a(1)',ia,ib,s2a,r5a,polq(2,1),abs(1.0 - c12*s2a - c22*s2a**2),&
      !abs(1.0 - c12*s2a - c22*s2a**2)*r5a*polq(2,1)
      do j = 1, 23
         aiw_a  = abs(1.0 - c11*s1a)               * r3a * pol(ia,j)
         aiw_b  = abs(1.0 - c11*s1b)               * r3b * pol(ib,j)
         aiw_out(j) = aiw_a + (aiw_b - aiw_a) *w3
         !aiw_out(j) = ((aiw_a/denom_3a) + (aiw_b/denom_3b))/((1/denom_3a)+(1/denom_3b))

         aiwq_a = abs(1.0 - c12*s2a - c22*s2a**2)  * r5a * polq(ia,j)
         aiwq_b = abs(1.0 - c12*s2b - c22*s2b**2)  * r5b * polq(ib,j)
         aiwq_out(j) = abs(aiwq_a + (aiwq_b - aiwq_a) *w5)
         !aiwq_out(j) = ((aiwq_a/denom_5a) + (aiwq_b/denom_5b))/((1/denom_5a)+(1/denom_5b))
      end do
      !write(5,*)'aiw_a,aiw_b',abs(1.0 - c11*s1a)*r3a*pol(ia,1),abs(1.0 - c11*s1b)*r3b*pol(ib,1)
      !write(5,*)'aiwq_a,aiwq_b',abs(1.0 - c12*s2a - c22*s2a**2)  * r5a *polq(ia,1),abs(1.0 - c12*s2b - c22*s2b**2)  &
      !* r5b * polq(ib,1)
      !write(5,*)'aiw_out(j),aiwq_out(j)',aiw_out(1),aiwq_out(1)

   end subroutine polar_weighted_avg

end module polar_interp
