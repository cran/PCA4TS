!======================================================================
!  This code is used to calculate the thresholding procedure for the segment procedure 
!======================================================================

!==============================================================
!                         INPUTS
!  SY: Sigma_y     Y: the observation matrix   my: mean_y
!  k: interger between 0 and k_0      n: the number of the observations 
!  p: the dimension of the time series   delta: the value of tuning parameter
!  b: the vector of the thresholding matrix
!==============================================================
subroutine thresh(SY,Y,my,k,n,p,delta,b)     
implicit double precision (a-h,q-z)
integer n,p,i,j,k,s,t
double precision a(p,p),Y(p,n),SY(p,p)
double precision my,theta,lam,p1,b,delta
dimension my(p),b(p*p)

s=0
p1=p
do i=1,p
    do j=1,p
    s=s+1
    theta=0
       do t=1,(n-k)
         theta=theta+((Y(i,t+k)-my(i))*(Y(j,t)-my(j))-SY(i,j))**2
       enddo
    theta=theta/n
    lam=delta*sqrt(theta*log(p1)/n)
      if(abs(SY(i,j)).lt.lam) then
      SY(i,j)=0
      endif
     b(s)=SY(i,j)
    enddo
enddo 



end
