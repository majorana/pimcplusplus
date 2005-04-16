      function p1 (t)
      implicit none
      real*8 p1,t
      p1 = ((t-1.0)*(t-1.0)*(1.0+2.0*t))
      end function p1

      function p2 (t)
      implicit none
      real*8 p2,t
      p2 = (t*t*(3.0-2.0*t))
      end function p2

      function q1 (t)
      implicit none
      real*8 q1,t
      q1 = (t*(t-1.0)*(t-1.0))
      end function q1

      function q2 (t)
      implicit none
      real*8 q2,t
      q2 = (t*t*(t-1.0))
      end function q2


      subroutine z3spline (x,y,z,x0,dx,nx,y0,dy,ny,z0,dz,nz,F,num,vals)
      implicit none
      real*8 x,x0,dx,y,y0,dy,z,z0,dz,xlo,ylo,zlo,u,v,w,
     +       a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3
      real*8 p1,p2,q1,q2
      integer nx,ny,nz, ixl,ixh,iyl,iyh,izl,izh,num,i
      complex*16 F(8,num,nz,ny,nx)
      complex*16 vals(num)
      
      ixl = int((x-x0)/dx)+1
      ixh = ixl+1
      xlo = dx*(ixl-1)
      iyl = int((y-y0)/dy)+1
      iyh = iyl+1
      ylo = dx*(iyl-1)
      izl = int((z-z0)/dx)+1
      izh = izl + 1
      zlo = dz*(izl-1)

      u = (x-xlo)/dx
      v = (y-ylo)/dy
      w = (z-zlo)/dz
      
      a0 = p1(u)
      a1 = p2(u)
      a2 = dx*q1(u)
      a3 = dx*q2(u)

      b0 = p1(v)
      b1 = p2(v)
      b2 = dy*q1(v)
      b3 = dy*q2(v)

      c0 = p1(w)
      c1 = p2(w)
      c2 = dz*q1(w)
      c3 = dz*q2(w)

      do i = 1,num
         vals(i) =a0*(b0*(c0*F(1,i,izl,iyl,ixl)+c1*F(1,i,izh,iyl,ixl)+
     +                    c2*F(4,i,izl,iyl,ixl)+c3*F(4,i,izh,iyl,ixl))+
     +                b1*(c0*F(1,i,izl,iyh,ixl)+c1*F(1,i,izh,iyl,ixl)+
     +                    c1*F(4,i,izl,iyh,ixl)+c3*F(4,i,izh,iyh,ixl))+
     +                b2*(c0*F(3,i,izl,iyl,ixl)+c1*F(3,i,izh,iyl,ixl)+
     +                    c2*F(7,i,izl,iyl,ixl)+c3*F(7,i,izh,iyl,ixl))+
     +                b3*(c0*F(3,i,izl,iyh,ixl)+c1*F(3,i,izh,iyl,ixl)+
     +                    c1*F(7,i,izl,iyh,ixl)+c3*F(7,i,izh,iyh,ixl)))+
     +            a1*(b0*(c0*F(1,i,izl,iyl,ixh)+c1*F(1,i,izh,iyl,ixh)+
     +                    c2*F(4,i,izl,iyl,ixh)+c3*F(4,i,izh,iyl,ixh))+
     +                b1*(c0*F(1,i,izl,iyh,ixh)+c1*F(1,i,izh,iyl,ixh)+
     +                    c1*F(4,i,izl,iyh,ixh)+c3*F(4,i,izh,iyh,ixh))+
     +                b2*(c0*F(3,i,izl,iyl,ixh)+c1*F(3,i,izh,iyl,ixh)+
     +                    c2*F(7,i,izl,iyl,ixh)+c3*F(7,i,izh,iyl,ixh))+
     +                b3*(c0*F(3,i,izl,iyh,ixh)+c1*F(3,i,izh,iyl,ixh)+
     +                    c1*F(7,i,izl,iyh,ixh)+c3*F(7,i,izh,iyh,ixh)))+
     +            a2*(b0*(c0*F(2,i,izl,iyl,ixl)+c1*F(2,i,izh,iyl,ixl)+
     +                    c2*F(6,i,izl,iyl,ixl)+c3*F(6,i,izh,iyl,ixl))+
     +                b1*(c0*F(2,i,izl,iyh,ixl)+c1*F(2,i,izh,iyl,ixl)+
     +                    c1*F(6,i,izl,iyh,ixl)+c3*F(6,i,izh,iyh,ixl))+
     +                b2*(c0*F(5,i,izl,iyl,ixl)+c1*F(5,i,izh,iyl,ixl)+
     +                    c2*F(8,i,izl,iyl,ixl)+c3*F(8,i,izh,iyl,ixl))+
     +                b3*(c0*F(5,i,izl,iyh,ixl)+c1*F(5,i,izh,iyl,ixl)+
     +                    c1*F(8,i,izl,iyh,ixl)+c3*F(8,i,izh,iyh,ixl)))+
     +            a3*(b0*(c0*F(2,i,izl,iyl,ixh)+c1*F(2,i,izh,iyl,ixh)+
     +                    c2*F(6,i,izl,iyl,ixh)+c3*F(6,i,izh,iyl,ixh))+
     +                b1*(c0*F(2,i,izl,iyh,ixh)+c1*F(2,i,izh,iyl,ixh)+
     +                    c1*F(6,i,izl,iyh,ixh)+c3*F(6,i,izh,iyh,ixh))+
     +                b2*(c0*F(5,i,izl,iyl,ixh)+c1*F(5,i,izh,iyl,ixh)+
     +                    c2*F(8,i,izl,iyl,ixh)+c3*F(8,i,izh,iyl,ixh))+
     +                b3*(c0*F(5,i,izl,iyh,ixh)+c1*F(5,i,izh,iyl,ixh)+
     +                    c1*F(8,i,izl,iyh,ixh)+c3*F(8,i,izh,iyh,ixh)))

      end do
     

      end subroutine z3spline



