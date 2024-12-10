## Fortran codes for height anomaly computation on or outside geoid using Stokes/Hotine numerical integral
https://www.zcyphygeodesy.com/en/h-nd-133.html
## [Algorithm purpose]
    Using the generalized Stokes/Hotine numerical integral, from the ellipsoidal height grid of the equipotential surface and gravity anomaly or disturbance (mGal) grid on the surface, compute the height anomaly (m) on or outside the geoid.
    Here the equipotential boundary surface does not need to be the geoid, which can be any equipotential (or normal/orthometric equiheight) surface. E.g., the ground residual height anomaly can be calculated from the residual gravity anomaly / disturbance on the equipotential surface outside the ground.
Height anomaly on the geoid is equal to the geoid undulation, that is, the geoidal (ellipsoidal) height.
    The Stokes boundary value theory requires that the boundary surface should be an equipotential surface, that is, the gravity anomaly/disturbance should be on the equipotential surface.
    It is usually necessary to employ the remove-restore scheme with a reference geopotential model to use the finite radius for gravity field integral. Firstly, remove model gravity anomaly/disturbance on the boundary surface, then integrate to obtain the residual height anomaly at the calculation point, and finally restore the model height anomaly at the calculation point.
    The equipotential surface can be constructed from a global geopotential model (not greater than 360 degrees), which can also be represent by a normal (orthometric) equiheight surface with the altitude of not more than ten kilometers.
## [Main program for test entrance]
    StokesHotinenumintegral.f90
    The record format of the input calculation point file: ID (point no / point name), longitude (decimal degrees), latitude (decimal degrees), ellipsoidal height (m)......
    Input the ellipsoidal height grid file of the equipotential boundary surface, which employed to calculate the integral distance.
    Input the residual gravity anomaly/disturbance (mGal) grid file on the equipotential surface.
    Input parameter mode - when mode =0 for Stokes integral, and when mode = 1 for Hotine integral.
    The record format of the output file reslt.txt: Behind the record of the calculation point file, appends 1 column of residual height anomaly (m).
## (1) Algorithm module for the generalized Stokes numerical integral
    Real*8 StokesBLH(BLH,gra,sfh,nlat,nlon,hd,dr,GRS)
    Input parameters: BLH(3) - longitude (decimal degrees), latitude (decimal degrees), ellipsoidal height (m) of the calculation point.
    Input parameters: sfh(nlat,nlon) - the ellipsoidal height grid of the equipotential boundary surface, which employed to calculate the integral distance.
    Input parameters: gra(nlat,nlon) - the residual gravity anomaly (mGal) grid on the equipotential surface.
    Input parameters: dr, hd(6) - the integral radius (m) and grid specification parameters (minimum and maximum longitude, minimum and maximum latitude, longitude and latitude intervals of a cell grid).
    Input parameters: GRS(6) - gm, ae, j2, omega, 1/f, default value
    Return: the calculated residual height anomaly(m).
## (2) Algorithm module for the generalized Stokes numerical integral
    Real*8 HotineBLH(BLH,rga,sfh,nlat,nlon,hd,dr,GRS)
    Input parameters: rga(nlat,nlon) - the residual gravity disturbance (mGal) grid on the equipotential surface.
    Return: the calculated residual height anomaly(m).
## (3) Calculation module for the normal gravity field
    normdjn(GRS,djn); GNormalfd(BLH,NFD,GRS)
    Return parameters: NFD(5) - the normal geopotential (m2/s2), normal gravity (mGal), normal gravity gradient (E), normal gravity line direction (', expressed by its north declination relative to the center of the Earth center of mass) or normal gravity gradient direction (', expressed by its north declination relative to the Earth center of mass).
## (4) Calculation module for Legendre functions and their derivatives to ψ
    LegPn_dt2(pn,dp1,dp2,n,t) ! t=cos ψ
## (5) Algorithm library for transforming of geodetic coordinates
    BLH_RLAT(GRS, BLH, RLAT); BLH_XYZ(GRS, BLH, XYZ)
    RLAT_BLH(GRS, RLAT, BLH)
## (6) Algorithm library for interpolation point value from numerical grid
    CGrdPntD(lon,lat,dt,row,col,hd); CGrdPntD2(lon,lat,dt,row,col,hd)
    CShepard(lon,lat,dt,row,col,hd); Gauss2D(lon,lat,dt,row,col,hd)
## (7) Other auxiliary modules
    PickRecord(str0, kln, rec, nn)
## [For compile and link]
    Fortran90, 132 Columns fixed format. Fortran compiler for any operating system. No external link library required.
## [Algorithmic formula] PAGravf4.5 User Reference https://www.zcyphygeodesy.com/en/
    1.4.1 Format convention for geodetic data file
    7.9.1 Stokes and Hotine integral formulas outside geoid
    7.1(4) Low-dgree Legendre function and its first and second derivative algorithms
DOS executable test file and all input and output data.
