# Main file for julia wrapper of NASA MarsGRAM 2010 (written in Fortran)
# Uses Interpolation.jl
# Relies on MarsGRAM2010 documentation

# OUTPUT : Obtain Mars Atmospheric density at any altitude value
#  atmospheric_density_interpolated_marsgram(altitude)


# TO BE MODIFIED ON EACH MACHINE
# Localize and open the .so file from marsgram
# using Libdl
# marsgram2 = Libdl.dlopen("/Users/kevintracy/devel/MarsGRAM2010/Code/lib.so")

const libso = "/Users/kevintracy/devel/MarsGRAM2010/Code/lib.so"
# Files containint Interpolation File
# include("interpolation.jl")

# Wrapper of MarsGram2010 main routine for Mars Density Profile.
# Default values of input variables are the same as the original code in Fortran.
# See meaning of Variables in MarsGRAM2010 documentation.

# Code below only calls setup_m10_ and datastep_m10_ routines from FORTRAN code

function wrapper_density(LSTFL="INPUT.txt", OUTFL="OUTPUT.txt",
        PROFILE="null", WAVEFILE="null", DATADIR ="null", GCMDIR="null",
        IERT=1,IUTC=1, MONTH=7, MDAY=20, MYEAR=20, NPOS=400, IHR=12, IMIN=30,
        SEC=0.0, LONEW=0, DUSTTAU=0, DUSTMIN=0.3, DUSTMAX=1.0, DUSTNU=0.003,
        DUSTDIAM=5.0, DUSTDENS=3000., ALS0=0.0, ALSDUR=48., INTENS=0.0,
        RADMAX=0.0, DUSTLAT=0.0, DUSTLON=0.0, F107=68.0, NR1=1234, NVARX=1, NVARY=0,
        LOGSCALE=0, FLAT=22.48, FLON=47.97, FHGT=0.0, MOLAHGTS=1, HGTASFCM=0.0,
        ZOFFSET=3.25, IBOUGHER=1, DELHGT=0.4, DELLAT=0.5, DELLON=0.5, DELTIME=500.0,
        DELTATEX=0.0, PROFNEAR=0.0, PROFFAR=0.0, RPSCALE=1.0, RWSCALE=1.0, WLSCALE=1.0,
        WMSCALE=1.0, BLWINFAC=1.0, NMONTE=1, IUP=13, WAVEA0=1.0, WAVEDATE=0.0, WAVEA1=0.0,
        WAVEPHI1=0.0, PHI1DOT=0.0, WAVEA2=0.0, WAVEPHI2=0.0, PHI2DOT=0.0, WAVEA3=0.0,
        WAVEPHI3=0.0, PHI3DOT=0.0, IUWAVE=0, WSCALE=20., CORLMIN=0.0, IPCLAT=1, REQUA=3396.19,
        RPOLE=3376.20, MAPYEAR=0, IDAYDATA=1)

        #Variable Initialization for setup subroutine
        CHGT = [1.0]; CLAT = [1.0]; CLON = [1.0]; CSEC = [1.0]; DAY0 = [1.0]; RHOD=[1.0];
        RHOU=[1.0]; RHOV=[1.0]; RHOW=[1.0]; DHGT=[1.0]; DLAT=[1.0]; DLON=[1.0]; DTIME=[1.0];
        MAXNUM = [1]; NRN1=[1]; NMCR1=[1]; DSUNLS=[1.0]; DRADAU=[1.0]; DOWLT=[1.0];
        LNEW=[1]; IUSTDOUT = 6; IULIST=[1]; HGTASFC=[1.0]; INERT=[1]; INUTC=[1]; STEPMIN=[1.0];
        PROFNR=[1.0]; PROFFR=[1.0]; NPROF=[1]

        #call setup fotran subroutine, setup variables for the run
        #Pass by references for Fortran code (match the variables types in fortran)
        # ccall(Libdl.dlsym(marsgram2,:setup_m10_),
        #     Nothing, (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
        #     Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
        #     Ref{Float64}, Ref{Float64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Float64},
        #     Ref{Float64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Int64}, Ref{Int64},
        #     Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64}, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8},
        #     Ref{Int64}, Ref{Int64},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Float64},Ref{Int64},Ref{Int64}, Ref{Float64},Ref{Float64},
        #     Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},
        #     Ref{Float64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64},
        #     Ref{Float64},Ref{Float64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},Ref{Float64}, Ref{Float64},
        #     Ref{Float64}, Ref{Float64},Ref{Float64}, Ref{Float64},Ref{Float64}, Ref{Float64},
        #     Ref{Int64}, Ref{Int64}, Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},
        #     Ref{Float64},Ref{Float64},Ref{Float64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Int64}, Ref{Float64}, Ref{Float64},
        #     Ref{Int64}, Ref{Int64}), CHGT, CLAT, CLON, CSEC, DAY0, RHOD, RHOU, RHOV, RHOW,
        #     DHGT, DLAT, DLON, DTIME, MAXNUM, NRN1, NMCR1, DSUNLS, DRADAU, DOWLT, LNEW, IUSTDOUT, IULIST, HGTASFC,
        #     INERT, INUTC, STEPMIN, PROFNR, PROFFR, NPROF, LSTFL, OUTFL, PROFILE, WAVEFILE, DATADIR, GCMDIR, IERT, IUTC, MONTH, MDAY, MYEAR, NPOS
        #     ,IHR, IMIN, SEC, LONEW, DUSTTAU, DUSTMIN, DUSTMAX, DUSTNU, DUSTDIAM, DUSTDENS
        #     ,ALS0, ALSDUR, INTENS, RADMAX, DUSTLAT, DUSTLON, F107, NR1, NVARX, NVARY, LOGSCALE
        #     ,FLAT, FLON, FHGT, MOLAHGTS, HGTASFCM, ZOFFSET, IBOUGHER, DELHGT, DELLAT, DELLON, DELTIME
        #     , DELTATEX, PROFNEAR, PROFFAR, RPSCALE, RWSCALE, WLSCALE, WMSCALE, BLWINFAC, NMONTE, IUP
        #     , WAVEA0, WAVEDATE, WAVEA1, WAVEPHI1, PHI1DOT, WAVEA2, WAVEPHI2, PHI2DOT, WAVEA3, WAVEPHI3
        #     ,PHI3DOT, IUWAVE, WSCALE, CORLMIN, IPCLAT, REQUA, RPOLE, MAPYEAR, IDAYDATA)
        ccall((:setup_m10_,libso),
            Nothing, (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
            Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
            Ref{Float64}, Ref{Float64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Float64},
            Ref{Float64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Int64}, Ref{Int64},
            Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64}, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8},
            Ref{Int64}, Ref{Int64},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Int64},Ref{Float64},Ref{Int64},Ref{Int64}, Ref{Float64},Ref{Float64},
            Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},
            Ref{Float64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64},
            Ref{Float64},Ref{Float64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},Ref{Float64}, Ref{Float64},
            Ref{Float64}, Ref{Float64},Ref{Float64}, Ref{Float64},Ref{Float64}, Ref{Float64},
            Ref{Int64}, Ref{Int64}, Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},
            Ref{Float64},Ref{Float64},Ref{Float64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Int64}, Ref{Float64}, Ref{Float64},
            Ref{Int64}, Ref{Int64}), CHGT, CLAT, CLON, CSEC, DAY0, RHOD, RHOU, RHOV, RHOW,
            DHGT, DLAT, DLON, DTIME, MAXNUM, NRN1, NMCR1, DSUNLS, DRADAU, DOWLT, LNEW, IUSTDOUT, IULIST, HGTASFC,
            INERT, INUTC, STEPMIN, PROFNR, PROFFR, NPROF, LSTFL, OUTFL, PROFILE, WAVEFILE, DATADIR, GCMDIR, IERT, IUTC, MONTH, MDAY, MYEAR, NPOS
            ,IHR, IMIN, SEC, LONEW, DUSTTAU, DUSTMIN, DUSTMAX, DUSTNU, DUSTDIAM, DUSTDENS
            ,ALS0, ALSDUR, INTENS, RADMAX, DUSTLAT, DUSTLON, F107, NR1, NVARX, NVARY, LOGSCALE
            ,FLAT, FLON, FHGT, MOLAHGTS, HGTASFCM, ZOFFSET, IBOUGHER, DELHGT, DELLAT, DELLON, DELTIME
            , DELTATEX, PROFNEAR, PROFFAR, RPSCALE, RWSCALE, WLSCALE, WMSCALE, BLWINFAC, NMONTE, IUP
            , WAVEA0, WAVEDATE, WAVEA1, WAVEPHI1, PHI1DOT, WAVEA2, WAVEPHI2, PHI2DOT, WAVEA3, WAVEPHI3
            ,PHI3DOT, IUWAVE, WSCALE, CORLMIN, IPCLAT, REQUA, RPOLE, MAPYEAR, IDAYDATA)

        # #initialize values
        # NUMWAVE = 0; PERTSTEP =[0.0]; IUPDATE = 0; EOF = 0; TEMP =[1.0]; PRES =[1.0];
        # DENSLO=[1.0];DENS=[1.0];DENSHI=[1.0]; DENSP=[1.0]; EWWIND=[1.0]; EWPERT=[1.0];
        # NSWIND=[1.0]; NSPERT=[1.0]; VWPERT=[1.0]; HRHO=[1.0]; HPRES=[1.0]; CORLIM=[1.0];
        # DENSTOT=[1.0]; ALS=[1.0]; SZA=[1.0]; OWLT=[1.0]; SUNLAT=[1.0]; SUNLON=[1.0];
        # MARSAU=[1.0]; TLOCAL=[1.0];
        #
        # # Initialize density array
        # dens_results = []
        # height = FHGT:DELHGT:FHGT+(NPOS-1)*DELHGT
        #
        # # Call datastep_m10_ fortran routine at each step. Returning values of
        # # atmospheric parameters at each step of the loop.
        # # Modifications expected if MonteCarlo Number is more than 1.
        # for I=0:1:MAXNUM[1]
        #     ccall(Libdl.dlsym(marsgram2,:datastep_m10_), Nothing, (Ref{Int64}, Ref{Float64},
        #     Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64},
        #     Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
        #     Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
        #     Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
        #     Ref{Float64}, Ref{Float64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Int64}, Ref{Float64}, Ref{Int64},
        #     Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
        #     Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64}), I, CHGT, CLAT, CLON, CSEC, DAY0,
        #     RHOD, RHOU, RHOV, RHOW, EOF, DELHGT, DELLAT, DELLON, DELTIME, TEMP, PRES,
        #     DENSLO, DENS, DENSHI, DENSP, EWWIND, EWPERT, NSWIND, NSPERT, VWPERT, HRHO, HPRES,
        #     0.0, 0.0, 0.0, 0.0, 0.0, LONEW, CORLIM, DENSTOT, NUMWAVE, HGTASFC, IERT, IUTC, PERTSTEP,
        #     CORLMIN, IUPDATE, ALS, SZA, OWLT, SUNLAT, SUNLON, MARSAU, TLOCAL, PROFNEAR, PROFFAR, NPROF)
        #     append!(dens_results, DENS)
        #     # Here only use the values of DENS (given by datastep_m10_ function)
        #     # Can be extended to use other information from MarsGRAM
        #     # e.g. PRES, DENSLOW, DENSHIGH etc...
        # end

        # return height, dens_results
        return nothing
end

# Test Functions

# Get the density profile from wrapper_density function (from fortran MarsGRAM)
# height, density = wrapper_density()
wrapper_density()

# Perform Chebyshev Polynomial Interpolation of MarsGram data
# interpol_coeff = compute_chebyshev_coefficients_atmosphere_log(height, density, 10)
#
#
# # Function below returns the interpolated Mars atmopsheric density value for any
# # value of altitude.
#
# function atmospheric_density_interpolated_marsgram(interpol_coeff, altitude)
#     #altitude in km
#     #density in kg.m-3
#     alt_min = minimum(height);
#     alt_max = maximum(height);
#     density = atmosphere_density_chebyshev(interpol_coeff, altitude, alt_min, alt_max)
#     return density
# end
