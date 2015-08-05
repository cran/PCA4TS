!======================================================================
! 
! This code is used to speed up the core part in segmenting the volatility process
! 
!======================================================================

subroutine volsegment(Y,meanY,n,p,k0,res)   
implicit double precision (a-h,q-z)
integer n,p,l,i,j,k,k0,t,s,m,i1,i2,i3
double precision Y(p,n),SigmaY(p,p),WY(p,p),SigmaY1(p,p)
double precision meanY,y2,y1,res1,res2,y3,res
dimension meanY(p),y2(p),res(p*p),y1(p),y3(p)

do i=1,p
   do j=1,p
      WY(i,j)=0
   enddo
enddo 

do l=1,n
  y2=Y(:,l)-meanY   
    do k=1,k0
       
       do i=1,p
         do j=1,p
           SigmaY(i,j)=0
          
         enddo
       enddo
      
       do t=(k+1),n
          s=t-k
          y1=Y(:,s)-meanY
          res1=0
          res2=0 
            do m=1,p
              res1=res1+y1(m)*y1(m)
              res2=res2+y2(m)*y2(m)
            enddo 
          if(res1<res2) then
            y3=Y(:,t)-meanY  
            do i1=1,p
              do i2=1,p
                SigmaY(i1,i2)=SigmaY(i1,i2)+y3(i1)*y3(i2)
              enddo
            enddo  
          endif
       enddo         
         Sigmay=Sigmay/(n-k)
         Sigmay=matmul(Sigmay,Sigmay)
         WY=WY+Sigmay
    enddo
enddo


s=0
do i=1,p
  do j=1,p
  s=s+1
  res(s)=WY(i,j)
  ! res(s)=Sigmay(i,j)
  enddo
enddo


end
