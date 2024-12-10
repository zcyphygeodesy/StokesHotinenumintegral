      real*8 function HotineBLH(BLH,gra,sfh,nlat,nlon,hd,dr,GRS)
      !the generalized Hotine rigorous numerical integral 广义Hotine积分积分
      !dr-积分半径 the integral radius (m)
!-------------------------------------------------------------
      implicit none
	integer::i,j,nlat,nlon,i0,j0,ni,nj
	real*8::dr,gra(nlat,nlon),sfh(nlat,nlon)
	real*8::hd(6),pi,RAD,ds,mdr,tt,rr,r0,r1,r2,tmp
	real*8::GRS(6),BLH(3),XYZ(3),rln(3),BLH0(3),XYZ0(3),BLH1(3),XYZ1(3),rln1(3)
	real*8 CGrdPntD2,HotineSgn,L0,L1,SK,rst,gr,NFD(5)
!-----------------------------------------------------------------
      pi=datan(1.d0)*4.d0;RAD=pi/180.d0;BLH0=BLH
      BLH0(3)=CGrdPntD2(BLH(2),BLH(1),sfh,nlat,nlon,hd)!计算点正下方地面点位置
      call BLH_RLAT(GRS,BLH,rln);rr=rln(1);call BLH_XYZ(GRS,BLH,XYZ)
      call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
      call BLH_XYZ(GRS,BLH0,XYZ0)
      r0=dsqrt(XYZ0(1)**2+XYZ0(2)**2+XYZ0(3)**2)
      mdr=r0*hd(5)*RAD*dcos(rln(2)*RAD)/2.d0 !奇异点判断
      ni=nint(dr/r0/RAD/hd(6)+1.d0) !积分半径dr对应的地面格网数
      nj=nint(dr/r0/RAD/hd(5)/dcos(rln(2)*RAD)+1.d0)
	i0=nint((BLH(1)-hd(3))/hd(6)+0.5d0)
	j0=nint((BLH(2)-hd(1))/hd(5)+0.5d0)!计算点所在的地面格网i0,j0
      rst=0.d0!gra0*(r0/rr)**2/gr*1.d-5!!特殊改进
	do i=i0-ni,i0+ni
	  if(i<1.or.i>nlat)goto 9100
        BLH1(1)=hd(3)+(real(i)-0.5d0)*hd(6)
	  do j=j0-nj,j0+nj
	    if(j<1.or.j>nlon)goto 9101
	    BLH1(2)=hd(1)+(real(j)-0.5d0)*hd(5)
          BLH1(3)=sfh(i,j);call BLH_XYZ(GRS,BLH1,XYZ1)
          L0=dsqrt((XYZ1(1)-XYZ0(1))**2+(XYZ1(2)-XYZ0(2))**2+(XYZ1(3)-XYZ0(3))**2)
          if(L0>dr)goto 9101
          L1=dsqrt((XYZ1(1)-XYZ(1))**2+(XYZ1(2)-XYZ(2))**2+(XYZ1(3)-XYZ(3))**2)
          if(L0<mdr)then!计算奇异积分
             tmp=HotineSgn(BLH,gra,sfh,nlat,nlon,hd,i,j,4,GRS)
             rst=rst+tmp
             goto 9101 
          endif
          call BLH_RLAT(GRS,BLH1,rln1);r1=rln1(1)
          tt=1.d0-2.d0*(L0/r1/2.d0)**2
          ds=hd(5)*hd(6)*RAD**2*dcos(rln1(2)*RAD)
          SK=2.d0/L1-1.d0/rr-3.d0*r1*tt/rr**2-dlog(dabs(rr+L1-r1*tt)/(rr*(1.d0-tt)))/r1
          rst=rst+r1**2*SK*gra(i,j)*ds/gr/4.d0/pi*1.d-5
9101      continue
	  enddo
9100    continue
	enddo
	HotineBLH=rst
9002	return
      end
!--------------------------------------------------------------------------------
      real*8 function HotineSgn(BLH,gra,sfh,nlat,nlon,hd,i0,j0,m,GRS)
      !细化核函数，计算BLH点的Hotine奇异积分
      !m-核函数细化为m*m
      !i0,j0-奇异点格网位置
!-------------------------------------------------------------
      implicit none
	integer::m,i,j,nlat,nlon,i0,j0,ni,nj
	real*8::dr,gra(nlat,nlon),sfh(nlat,nlon)
	real*8::hd(6),pi,RAD,ds,mdr,tt,rr,r0,r1,r2,rst,rv
	real*8::GRS(6),BLH(3),XYZ(3),rln(3),BLH0(3),XYZ0(3),rln0(3),BLH1(3),XYZ1(3),rln1(3)
	real*8 CGrdPntD2,L0,L1,SK,lon,lat,dg,gr,NFD(5),PS
!-----------------------------------------------------------------
      BLH0=BLH;rv=hd(5)/dble(m);pi=datan(1.d0)*4.d0;RAD=pi/180.d0
      lat=hd(3)+real(i0-1)*hd(6);lon=hd(1)+real(j0-1)*hd(5)!格网左下角经纬度
      BLH0(3)=CGrdPntD2(BLH(2),BLH(1),sfh,nlat,nlon,hd)!计算点正下方地面点位置
      call BLH_XYZ(GRS,BLH,XYZ);call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
      rr=dsqrt(XYZ(1)**2+XYZ(2)**2+XYZ(3)**2)
      call BLH_XYZ(GRS,BLH0,XYZ0);call BLH_RLAT(GRS,BLH0,rln0);r0=rln0(1)
      rst=0.d0;mdr=r0*rv*RAD*dcos(BLH(1)*RAD)/dble(m)/4.d0  !奇异点判断
	do i=1,m
        BLH1(1)=lat+(real(i)-0.5d0)*rv
	  do j=1,m
	    BLH1(2)=lon+(real(j)-0.5d0)*rv
          BLH1(3)=CGrdPntD2(BLH1(2),BLH1(1),sfh,nlat,nlon,hd)
          dg=CGrdPntD2(BLH1(2),BLH1(1),gra,nlat,nlon,hd)
          call BLH_XYZ(GRS,BLH1,XYZ1)
          L0=dsqrt((XYZ1(1)-XYZ0(1))**2+(XYZ1(2)-XYZ0(2))**2+(XYZ1(3)-XYZ0(3))**2)
          L1=dsqrt((XYZ1(1)-XYZ(1))**2+(XYZ1(2)-XYZ(2))**2+(XYZ1(3)-XYZ(3))**2)
          call BLH_RLAT(GRS,BLH1,rln1);r1=rln1(1)
          ds=rv**2*RAD**2*dcos(rln1(2)*RAD)
          if(L0<mdr)L0=mdr
          tt=dsqrt(1.d0-(L0/r1)**2)
          SK=2.d0/L1-1.d0/rr-3.d0*r1*tt/rr**2-dlog(dabs(rr+L1-r1*tt)/(rr*(1.d0-tt)))/r1
          rst=rst+r1**2*SK*dg*ds/gr/4.d0/pi*1.d-5
	  enddo
9100    continue
	enddo
	HotineSgn=rst
9002	return
      end
