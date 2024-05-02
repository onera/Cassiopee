! ================================================================
  module sort2M
! ================================================================

  implicit none

  private
  public :: sort2

  contains


! ===============================================================
  subroutine sort2(n,arr)
! ===============================================================
! sorts an array arr(2,1:n) into ascending order using Quicksort
! using arr(1,:) as the primary key and arr(2,:) as the secondary
! key

  implicit none

  integer,parameter         :: m=7,nstack=10000
  integer,intent(in)        :: n
  integer,dimension(2,n)    :: arr
  integer,dimension(nstack) :: istack
  integer                   :: i,ir,j,jstack,k,l
  integer                   :: a1,a2,temp


  continue


  jstack=0
  l=1
  ir=n

 1 if (ir-l.lt.M) then

     do j=l+1,ir
       a1=arr(1,j)
       a2=arr(2,j)
       do i=j-1,1,-1
         if (arr(1,i).lt.a1) goto 2
         if ((arr(1,i).eq.a1).and.(arr(2,i).lt.a2)) goto 2
         arr(1,i+1)=arr(1,i)
         arr(2,i+1)=arr(2,i)
       end do
       i=0
  2    arr(1,i+1)=a1
       arr(2,i+1)=a2
     end do

     if (jstack == 0) return
     ir = istack(jstack)
     l = istack(jstack-1)
     jstack = jstack-2

  else

     k=(l+ir)/2

     temp=arr(1,k)
     arr(1,k)=arr(1,l+1)
     arr(1,l+1)=temp

     temp=arr(2,k)
     arr(2,k)=arr(2,l+1)
     arr(2,l+1)=temp

     if ((arr(1,l+1).gt.arr(1,ir)).or. & 
         ((arr(1,l+1).eq.arr(1,ir)).and.(arr(2,l+1).gt.arr(2,ir)))) then
       temp=arr(1,l+1)
       arr(1,l+1)=arr(1,ir)
       arr(1,ir)=temp
       temp=arr(2,l+1)
       arr(2,l+1)=arr(2,ir)
       arr(2,ir)=temp
     end if

     if ((arr(1,l).gt.arr(1,ir)).or. & 
         ((arr(1,l).eq.arr(1,ir)).and.(arr(2,l).gt.arr(2,ir)))) then
       temp=arr(1,l)
       arr(1,l)=arr(1,ir)
       arr(1,ir)=temp
       temp=arr(2,l)
       arr(2,l)=arr(2,ir)
       arr(2,ir)=temp
     end if

     if ((arr(1,l+1).gt.arr(1,l)).or. & 
         ((arr(1,l+1).eq.arr(1,l)).and.(arr(2,l+1).gt.arr(2,l)))) then
       temp=arr(1,l+1)
       arr(1,l+1)=arr(1,l)
       arr(1,l)=temp
       temp=arr(2,l+1)
       arr(2,l+1)=arr(2,l)
       arr(2,l)=temp
     end if

     i=l+1
     j=ir
     a1=arr(1,l)
     a2=arr(2,l)
  3  continue
     i=i+1
     if (arr(1,i).lt.a1) goto 3
     if ((arr(1,i).eq.a1).and.(arr(2,i).lt.a2)) goto 3
  4  continue
     j=j-1
     if (arr(1,j).gt.a1) goto 4
     if ((arr(1,j).eq.a1).and.(arr(2,j).gt.a2)) goto 4
     if (j.lt.i) goto 5

     temp=arr(1,i)
     arr(1,i)=arr(1,j)
     arr(1,j)=temp
     temp=arr(2,i)
     arr(2,i)=arr(2,j)
     arr(2,j)=temp

     goto 3
  5  arr(1,l)=arr(1,j)
     arr(2,l)=arr(2,j)
     arr(1,j)=a1
     arr(2,j)=a2

     jstack=jstack+2
!    push pointers to larger subarray on stack, process smaller 
!    subarray immediately.
     if (jstack.gt.NSTACK) stop 'NSTACK too small in sort2'
     if (ir-1+1.ge.j-1)then
       istack(jstack)=ir
       istack(jstack-1)=i
       ir=j-1
     else
       istack(jstack)=j-1
       istack(jstack-1)=l
       l=i
    end if
  end if
  goto 1

  end subroutine sort2

  end module sort2M
