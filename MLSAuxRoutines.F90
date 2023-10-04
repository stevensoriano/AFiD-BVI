!-------------------------------------------------------
!     contains auxiallary routine for MLS interpolation
!-------------------------------------------------------

      subroutine compinvA(pinvA,invA)
      IMPLICIT NONE

      double precision, intent(in) ::  pinvA(4,4)
      double precision, intent(out) :: invA(4,4)
      double precision det, usdet

!VS     matrix inversion of 4x4 matrix taken from graphics
!VS     library of MESA3D 

      invA(1,1)= pinvA(2,2)*pinvA(3,3)*pinvA(4,4) - &
                 pinvA(2,2)*pinvA(4,3)*pinvA(3,4) - &
                 pinvA(2,3)*pinvA(3,2)*pinvA(4,4) + &
                 pinvA(2,3)*pinvA(4,2)*pinvA(3,4) + &
                 pinvA(2,4)*pinvA(3,2)*pinvA(4,3) - &
                 pinvA(2,4)*pinvA(4,2)*pinvA(3,3)
      
      invA(1,2)=-pinvA(1,2)*pinvA(3,3)*pinvA(4,4) + &
                 pinvA(1,2)*pinvA(4,3)*pinvA(3,4) + &
                 pinvA(1,3)*pinvA(3,2)*pinvA(4,4) - &
                 pinvA(1,3)*pinvA(4,2)*pinvA(3,4) - &
                 pinvA(1,4)*pinvA(3,2)*pinvA(4,3) + &
                 pinvA(1,4)*pinvA(4,2)*pinvA(3,3)
      
      invA(1,3)= pinvA(1,2)*pinvA(2,3)*pinvA(4,4) - &
                 pinvA(1,2)*pinvA(4,3)*pinvA(2,4) - &
                 pinvA(1,3)*pinvA(2,2)*pinvA(4,4) + &
                 pinvA(1,3)*pinvA(4,2)*pinvA(2,4) + &
                 pinvA(1,4)*pinvA(2,2)*pinvA(4,3) - &
                 pinvA(1,4)*pinvA(4,2)*pinvA(2,3)
      
      invA(1,4)=-pinvA(1,2)*pinvA(2,3)*pinvA(3,4) + &
                 pinvA(1,2)*pinvA(3,3)*pinvA(2,4) + &
                 pinvA(1,3)*pinvA(2,2)*pinvA(3,4) - &
                 pinvA(1,3)*pinvA(3,2)*pinvA(2,4) - &
                 pinvA(1,4)*pinvA(2,2)*pinvA(3,3) + &
                 pinvA(1,4)*pinvA(3,2)*pinvA(2,3)
 
      invA(2,1)=-pinvA(2,1)*pinvA(3,3)*pinvA(4,4) + &
                 pinvA(2,1)*pinvA(4,3)*pinvA(3,4) + &
                 pinvA(2,3)*pinvA(3,1)*pinvA(4,4) - &
                 pinvA(2,3)*pinvA(4,1)*pinvA(3,4) - &
                 pinvA(2,4)*pinvA(3,1)*pinvA(4,3) + &
                 pinvA(2,4)*pinvA(4,1)*pinvA(3,3)
 
      invA(2,2)= pinvA(1,1)*pinvA(3,3)*pinvA(4,4) - &
                 pinvA(1,1)*pinvA(4,3)*pinvA(3,4) - &
                 pinvA(1,3)*pinvA(3,1)*pinvA(4,4) + &
                 pinvA(1,3)*pinvA(4,1)*pinvA(3,4) + &
                 pinvA(1,4)*pinvA(3,1)*pinvA(4,3) - &
                 pinvA(1,4)*pinvA(4,1)*pinvA(3,3)
 
      invA(2,3)=-pinvA(1,1)*pinvA(2,3)*pinvA(4,4) + &
                 pinvA(1,1)*pinvA(4,3)*pinvA(2,4) + &
                 pinvA(1,3)*pinvA(2,1)*pinvA(4,4) - &
                 pinvA(1,3)*pinvA(4,1)*pinvA(2,4) - &
                 pinvA(1,4)*pinvA(2,1)*pinvA(4,3) + &
                 pinvA(1,4)*pinvA(4,1)*pinvA(2,3)
 
      invA(2,4)= pinvA(1,1)*pinvA(2,3)*pinvA(3,4) - &
                 pinvA(1,1)*pinvA(3,3)*pinvA(2,4) - &
                 pinvA(1,3)*pinvA(2,1)*pinvA(3,4) + &
                 pinvA(1,3)*pinvA(3,1)*pinvA(2,4) + &
                 pinvA(1,4)*pinvA(2,1)*pinvA(3,3) - &
                 pinvA(1,4)*pinvA(3,1)*pinvA(2,3)

      invA(3,1)= pinvA(2,1)*pinvA(3,2)*pinvA(4,4) - &
                 pinvA(2,1)*pinvA(4,2)*pinvA(3,4) - &
                 pinvA(2,2)*pinvA(3,1)*pinvA(4,4) + &
                 pinvA(2,2)*pinvA(4,1)*pinvA(3,4) + &
                 pinvA(2,4)*pinvA(3,1)*pinvA(4,2) - &
                 pinvA(2,4)*pinvA(4,1)*pinvA(3,2)

      invA(3,2)=-pinvA(1,1)*pinvA(3,2)*pinvA(4,4) + &
                 pinvA(1,1)*pinvA(4,2)*pinvA(3,4) + &
                 pinvA(1,2)*pinvA(3,1)*pinvA(4,4) - &
                 pinvA(1,2)*pinvA(4,1)*pinvA(3,4) - &
                 pinvA(1,4)*pinvA(3,1)*pinvA(4,2) + &
                 pinvA(1,4)*pinvA(4,1)*pinvA(3,2)

      invA(3,3)= pinvA(1,1)*pinvA(2,2)*pinvA(4,4) - &
                 pinvA(1,1)*pinvA(4,2)*pinvA(2,4) - &
                 pinvA(1,2)*pinvA(2,1)*pinvA(4,4) + &
                 pinvA(1,2)*pinvA(4,1)*pinvA(2,4) + &
                 pinvA(1,4)*pinvA(2,1)*pinvA(4,2) - &
                 pinvA(1,4)*pinvA(4,1)*pinvA(2,2)

      invA(3,4)=-pinvA(1,1)*pinvA(2,2)*pinvA(3,4) + &
                 pinvA(1,1)*pinvA(3,2)*pinvA(2,4) + &
                 pinvA(1,2)*pinvA(2,1)*pinvA(3,4) - &
                 pinvA(1,2)*pinvA(3,1)*pinvA(2,4) - &
                 pinvA(1,4)*pinvA(2,1)*pinvA(3,2) + &
                 pinvA(1,4)*pinvA(3,1)*pinvA(2,2)

      invA(4,1)=-pinvA(2,1)*pinvA(3,2)*pinvA(4,3) + &
                 pinvA(2,1)*pinvA(4,2)*pinvA(3,3) + &
                 pinvA(2,2)*pinvA(3,1)*pinvA(4,3) - &
                 pinvA(2,2)*pinvA(4,1)*pinvA(3,3) - &
                 pinvA(2,3)*pinvA(3,1)*pinvA(4,2) + &
                 pinvA(2,3)*pinvA(4,1)*pinvA(3,2)

      invA(4,2)= pinvA(1,1)*pinvA(3,2)*pinvA(4,3) - &
                 pinvA(1,1)*pinvA(4,2)*pinvA(3,3) - &
                 pinvA(1,2)*pinvA(3,1)*pinvA(4,3) + &
                 pinvA(1,2)*pinvA(4,1)*pinvA(3,3) + &
                 pinvA(1,3)*pinvA(3,1)*pinvA(4,2) - &
                 pinvA(1,3)*pinvA(4,1)*pinvA(3,2)

      invA(4,3)=-pinvA(1,1)*pinvA(2,2)*pinvA(4,3) + &
                 pinvA(1,1)*pinvA(4,2)*pinvA(2,3) + &
                 pinvA(1,2)*pinvA(2,1)*pinvA(4,3) - &
                 pinvA(1,2)*pinvA(4,1)*pinvA(2,3) - &
                 pinvA(1,3)*pinvA(2,1)*pinvA(4,2) + &
                 pinvA(1,3)*pinvA(4,1)*pinvA(2,2)
 
      invA(4,4)= pinvA(1,1)*pinvA(2,2)*pinvA(3,3) - &
                 pinvA(1,1)*pinvA(3,2)*pinvA(2,3) - &
                 pinvA(1,2)*pinvA(2,1)*pinvA(3,3) + &
                 pinvA(1,2)*pinvA(3,1)*pinvA(2,3) + &
                 pinvA(1,3)*pinvA(2,1)*pinvA(3,2) - &
                 pinvA(1,3)*pinvA(3,1)*pinvA(2,2)

      det = pinvA(1,1)*invA(1,1) + &
            pinvA(2,1)*invA(1,2) + &
            pinvA(3,1)*invA(1,3) + &
            pinvA(4,1)*invA(1,4)

      usdet = 1.0d0/det
      invA = invA/det
      
      return
      end
      
!-------------------------------------------------------
      subroutine invert4(a,ia)
       implicit none
       real :: a (4,4)
       real :: ia (4,4)
       real :: det_l,invdet
      
       det_l =  a(1,1)*(a(2,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))   &
                  -a(2,3)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))   &
                  +a(2,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2)))  &
          -a(1,2)*(a(2,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))   &
                  -a(2,3)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))   &
                  +a(2,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1)))  &
          +a(1,3)*(a(2,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))   &
                  -a(2,2)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))   &
                  +a(2,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))  &
          -a(1,4)*(a(2,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))   &
                  -a(2,2)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))   &
                  +a(2,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))

       ia(1,1) =  a(2,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))-a(2,3)* &
      (a(3,2)*a(4,4)-a(3,4)*a(4,2))+a(2,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))

       ia(2,1) = -a(2,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))+a(2,3)* &
      (a(3,1)*a(4,4)-a(3,4)*a(4,1))-a(2,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))

       ia(3,1) =  a(2,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))-a(2,2)* &
      (a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(2,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))

       ia(4,1) = -a(2,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))+a(2,2)* &
      (a(3,1)*a(4,3)-a(3,3)*a(4,1))-a(2,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))


       ia(1,2) = -a(1,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))+a(1,3)* &
      (a(3,2)*a(4,4)-a(3,4)*a(4,2))-a(1,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))

       ia(2,2) =  a(1,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))-a(1,3)* &
      (a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(1,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))

       ia(3,2) = -a(1,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))+a(1,2)* &
      (a(3,1)*a(4,4)-a(3,4)*a(4,1))-a(1,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))

       ia(4,2) =  a(1,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))-a(1,2)* &
      (a(3,1)*a(4,3)-a(3,3)*a(4,1))+a(1,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))


       ia(1,3) =  a(1,2)*(a(2,3)*a(4,4)-a(2,4)*a(4,3))-a(1,3)* &
      (a(2,2)*a(4,4)-a(2,4)*a(4,2))+a(1,4)*(a(2,2)*a(4,3)-a(2,3)*a(4,2))

       ia(2,3) = -a(1,1)*(a(2,3)*a(4,4)-a(2,4)*a(4,3))+a(1,3)* &
      (a(2,1)*a(4,4)-a(2,4)*a(4,1))-a(1,4)*(a(2,1)*a(4,3)-a(2,3)*a(4,1))

       ia(3,3) =  a(1,1)*(a(2,2)*a(4,4)-a(2,4)*a(4,2))-a(1,2)* &
      (a(2,1)*a(4,4)-a(2,4)*a(4,1))+a(1,4)*(a(2,1)*a(4,2)-a(2,2)*a(4,1))

       ia(4,3) = -a(1,1)*(a(2,2)*a(4,3)-a(2,3)*a(4,2))+a(1,2)* &
      (a(2,1)*a(4,3)-a(2,3)*a(4,1))-a(1,3)*(a(2,1)*a(4,2)-a(2,2)*a(4,1))


       ia(1,4) = -a(1,2)*(a(2,3)*a(3,4)-a(2,4)*a(3,3))+a(1,3)* &
      (a(2,2)*a(3,4)-a(2,4)*a(3,2))-a(1,4)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))

       ia(2,4) =  a(1,1)*(a(2,3)*a(3,4)-a(2,4)*a(3,3))-a(1,3)* &
      (a(2,1)*a(3,4)-a(2,4)*a(3,1))+a(1,4)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))

       ia(3,4) = -a(1,1)*(a(2,2)*a(3,4)-a(2,4)*a(3,2))+a(1,2)* &
      (a(2,1)*a(3,4)-a(2,4)*a(3,1))-a(1,4)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))

       ia(4,4) =  a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))-a(1,2)* &
      (a(2,1)*a(3,3)-a(2,3)*a(3,1))+a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))

      invdet = 1.0d0/det_l

      ia(:,:)=ia(:,:)*invdet


      return
      end

!-------------------------------------------------------

      subroutine inverseLU(a,c)
      implicit none
      real :: a(4,4), c(4,4)
      real :: L(4,4), U(4,4)
      real :: b(4), d(4), x(4)
      real :: coeff
      integer i,j,k

      L = 0.0 ; U = 0.0 ; b = 0.0
      
      !forward elimination
      do k=1,3
       do i=k+1,4
            coeff=a(i,k)/a(k,k)
            L(i,k) = coeff
            do j=k+1,4
                  a(i,j)=a(i,j)-coeff*a(k,j)
            end do
       end do
      end do

      !prepare L U
      do i=1,4
       L(i,i) = 1.0
      end do
      do j=1,4
       do i=1,j
            U(i,j) = a(i,j)
       end do
      end do

      !compute columns of c
      do k=1,4
       b(k)=1.0
       d(1)=b(1)
        do i=2,4
        d(i)=b(i)
         do j=1,i-1
            d(i) = d(i)-L(i,j)*d(j)
         end do
        end do
       !solve ux=d with back subs.
       x(4)=d(4)/U(4,4)
       do i=3,1,-1
        x(i) = d(i)
        do j=4,i+1,-1
          x(i)=x(i)-U(i,j)*x(j)
        end do
        x(i) = x(i)/u(i,i)
       end do
      !fill solns of x(n) to k of C
      do i=1,4
            c(i,k)=x(i)
      end do
            b(k) = 0.0
      end do

      return
      end
