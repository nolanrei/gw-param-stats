###############################################################################################
## Calculates z correlations on the tropical WACCM dataset
###############################################################################################
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sg

dat1 = nc.Dataset("/glade/derecho/scratch/nreilly/Input1_year1.nc","r")
dat2 = nc.Dataset("/glade/derecho/scratch/nreilly/Input2_year1.nc","r")
dat3 = nc.Dataset("/glade/derecho/scratch/nreilly/Input3_year1.nc","r")
dat4 = nc.Dataset("/glade/derecho/scratch/nreilly/Input4_year1.nc","r")

nlon = dat1.dimensions["lon"].size
nlat = dat1.dimensions["lat"].size
nlev = dat1.dimensions["lev"].size
ntim = dat1.dimensions["time"].size

s = 24      # stride
nchunks = nlon//s
ugw_col_corr = np.zeros((nlev, nlev))
# 85:106 is the -10:10 latitude band
for j in range(nchunks):
    #          orographic                                          frontal                                              convective
    ugw = dat3.variables["UTGWORO"][:,:,85:106,j*s:(j+1)*s] + dat3.variables["UTGWSPEC"][:,:,85:106,j*s:(j+1)*s] + dat4.variables["BUTGWSPEC"][:,:,85:106,j*s:(j+1)*s]
    # adding dimensions so the shapes are (ntim, nlev, 1, nlat, s)x(ntim, 1, nlev, nlat, s) and then average all but (nlev, nlev)
    ugw = np.expand_dims(ugw,axis=2)
    tmp = np.mean(ugw*ugw.reshape(ntim,1,nlev,21,s),axis=(0,3,4))
    ugw_col_corr += tmp/nchunks

# not normalized here -- I divide in the plotting code by sqrt(column norm)*sqrt(row norm)

zcorr = nc.Dataset("waccm-zcorr-tropical.nc","w")
z = zcorr.createDimension("z",nlev)
zp = zcorr.createDimension("zp",nlev)
ugw_corr = zcorr.createVariable("zcorr","f8",("z","zp"))
ugw_corr[:] = ugw_col_corr
zcorr.close()