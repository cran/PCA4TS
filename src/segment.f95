!====================================================
!  This code is used to calculate the Sigma_y  
!====================================================

!==============================================
!                      INPUTS
!         Y: the observations(p by n);  my=mean_y 
!         k: interger between 0 and k0
!         n: number of observations
!         p: the dimension of y_t
!                      OUTPUTS
!         b: a vector of the elements of Sigma_y
!==============================================

subroutine segment(Y,my,k,n,p,b)    
implicit double precision (a-h,q-z)
integer n,p,i,j,k,t,s
double precision a(p,p),Y(p,n)
double precision my,C,D,D1,b
dimension my(p),C(p),D(p),b(p*p)

do t=1,(n-k)
  s=0
  C=Y(:,(t+k))-my
  D=Y(:,t)-my
   do i=1,p
    do j=1,p
      s=s+1
      b(s)=b(s)+C(i)*D(j)  
     enddo
   enddo 
enddo
b=b/n
end
