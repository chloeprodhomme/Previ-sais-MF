#!/usr/bin/env python
# Read data from an opendap server
import netCDF4
from cdo import *
#very important to have this line of code to avoir to fill /tmp
cdo = Cdo(tempdir='/home/prodhommec/tmp/') 
import requests
import numpy as np
import numpy.ma as ma
from datetime import date, datetime, timedelta
from dateutil.relativedelta import relativedelta
import math
from glob import glob
from netCDF4 import num2date, date2num
from  xarray import DataArray

def lon_index(longitude, lon_bnds):
    """
    Fonction to find the index of longitude corresponding to lon_bnds
    issue: it does not handle properly file where the separation is not done at greenwitch
    longitude: np array of longitude
    lon_bnds: longitude boundaries for the box to select
    return: tuple containing 1 or 2 numpy array of index of the longitude within the box
            tuple of size 2, if the longitude are over greenwitch, size 1 otherwise            
    """
    #if longitude are in -180 move to 0 - 360
    lon_bnds=np.array(lon_bnds)
    lons=np.array(longitude, copy=True)
    #print lons
    #print True in list(lons<0)
    if True in list(lons<0):
        lons[lons<0]=lons[lons<0]+360
    #print lons
    #if indices in -180, 180 move to 0 -360
    print(lon_bnds != "all", lon_bnds)
    if lon_bnds != "all":
        if True in list(lon_bnds<0):
            lon_bnds[lon_bnds<0]=lon_bnds[lon_bnds<0]+360
    
    #if 
    
    #check if box is over separation (most of the time greenwitch pero sometimes longitude can be splitted somewhere else)
    #not done.... Problem for grid_T
    #poorly done do not know how to handle it
        if lon_bnds[0]>lon_bnds[1]:
        #return a list with index before greenwitch and indexes after
        
            lon_inds = [np.where((lons >= lon_bnds[0]))[0] , np.where((lons <= lon_bnds[1]))[0]]
            if list(lon_inds[1])==[]:
                lon_inds=[lon_inds[0]]
            if list(lon_inds[0])==[]:
                lon_inds=[lon_inds[1]]
        else:
            lon_inds = [np.where((lons >= lon_bnds[0]) & (lons <= lon_bnds[1]))[0]]
    else:
        lon_inds=[np.arange(lons.shape[0])]
    return(lon_inds)

def lonlat_index(latitude, longitude, lat_bnds, lon_bnds):
    """
    Fonction to find the index of longitude and latitude corresponding to lon_bnds and lat_bnds
    can handle 2D arrays, but the selection of longitude will be done based on the center of the latitude
    issue: it does not handle properly file where the separation is not done at greenwitch (grid_T)
    latitude: np array of latitude
    longitude: np array of longitude
    lon_bnds: longitude boundaries for the box to select
    lat_bnds: latitude boundaries for the box to select
    return: numpy array of longitude, latitude, corresponding latitude index,    corresponding longitude index         
    """
    #handle 2D latitude array
    
    if len(latitude.shape)==2:
        lat1D=np.array(latitude[:,0], copy=True)
    else:
        lat1D=np.array(latitude, copy=True)
    
    #print lats
    #print lat_bnds
    #print np.where((lats >= lat_bnds[0]) & (lats <= lat_bnds[1]))[0]
    lat_inds = np.where((lat1D >= lat_bnds[0]) & (lat1D <= lat_bnds[1]))[0]
    
    #handle 2D longitude array (we base the longitude selection on the center of the latitude box)
    if len(longitude.shape)==2:
        centerlat=lat_inds[len(lat_inds)//2]
        #print(lat_inds[len(lat_inds)/2])
        lon1D=np.array(longitude[centerlat,:], copy=True)
        #print lons
    else:
        lon1D=np.array(longitude, copy=True)
        
    lon_inds = lon_index(lon1D, lon_bnds)  
    #print lon_inds 
    
    return(lat1D, lon1D, lat_inds, lon_inds)    


def area_grid(lat, lon):
    """
    Calculate the area of each grid cell
    Area is in square meters

    Input
    -----------
    lat: vector of latitude in degrees
    lon: vector of longitude in degrees

    Output
    -----------
    area: grid-cell area in square-meters with dimensions, [lat,lon]

    Notes
    -----------
    adapted from on the function in
    https://github.com/chadagreene/CDT/blob/master/cdt/cdtarea.m
    """

    lon[lon>=180.0]=lon[lon>=180.0]-360

    xlon, ylat = np.meshgrid(lon, lat)
    #R = earth_radius(ylat)
    R=6378137
    dlat = abs(np.deg2rad(np.gradient(ylat, axis=0)))
    dlon = abs(np.deg2rad(np.gradient(xlon%360, axis=1)))
    #print(dlon)

    if np.any(dlon>(300*3.14/180)):
        print("issue with jumps in longitude")
        sys.exit(1)

    dy = dlat * R
    dx = dlon * R * np.cos(np.deg2rad(ylat))
    #print(dx)
    #print(dy)
    #print(dx)

    area = dy * dx

    #xda = DataArray(
    #    area,
    #    dims=["latitude", "longitude"],
    #    coords={"latitude": lat, "longitude": lon},
    #    attrs={
    #        "long_name": "area_per_pixel",
    #        "description": "area per pixel",
    #        "units": "m^2",
    #    },
    #)
    return area
#obsreg=obs

def area_av(array, pos_lat, pos_lon, lats, lons, opt="mean", weightcalc=True, weights2D=None):
    """
    make the regional average of an array accounting for the size of the grid point (evolving with latitude)
    array: np array with at least 2dimensions for latitude and longitude 
    pos_lat: position of the dimension latitude in array
    pos_lon: position of the dimension longitude in array
    lats: latitude corresponding to array
    lons: longitude corresponding to array
    opt: "mean" copute the regional mean /"sum": compute the weighted sum
    weightcalc: if True compute the weights using the fonction area_grid, else the weight should be provided with weight2D argument
    weights2D: 2D numpy array of weight correspond to the grid
    """
    
    dim=array.shape
    if weightcalc:
        weights2D = area_grid(lats, lons)
    #weights=np.swapaxes(extend_table(area_grid(lats, lons), np.delete(np.delete(dim, pos_lon), pos_lat), len(dim)-1, pos_lat))
    
    weights=extend_table(weights2D, np.delete(np.delete(dim, pos_lon), pos_lat))
    #plt.imshow(area_grid(lats, lons))
    weights=np.ma.array(weights, mask=array.mask)
    #print(area_grid(lats, lons).shape)
    #print(array.shape)
    #print(weights.shape)

    if opt=="mean":
        sumweigth=np.ma.sum(np.ma.sum(weights, axis=pos_lon),axis=pos_lat)
        array_av = np.ma.sum(np.ma.sum(weights*array, axis=pos_lon),axis=pos_lat)/sumweigth
    if opt=="sum":
        array_av = np.ma.sum(np.ma.sum(weights*array, axis=pos_lon),axis=pos_lat)

    return(array_av)
    
def createdatelst(sdate1, sdate2, smonlstint):
    sdatelst=[]
    sdate=sdate1
    while sdate<sdate2:
        sdatelst.append([sdate+relativedelta(months=+(mon-1)) for mon in smonlstint])
        sdate=sdate+relativedelta(months=+12)
    sdatelst=np.ndarray.flatten(np.array(sdatelst))
    return(sdatelst)
  

def extend_table(array, dims_expend):
    array_ext=array
    dims_expend=list(dims_expend)
    dims_expend.reverse()
    for dim in dims_expend:
        array_ext=np.expand_dims(array_ext, axis=0).repeat(dim, axis=0)
        
    return(array_ext)


def getvarlat(varf):
    for varlat in ["Y", "lat", "latitude", "nav_lat"]:
        if varlat in varf.variables.keys():
            return varlat
        
def getvarlon(varf):
    for varlat in ["X", "lon", "longitude", "nav_lon"]:
        if varlat in varf.variables.keys():
            return varlat

def getvarmask(varf):
    for varmask in ["LSM", "land"]:
        if varmask in varf.variables.keys():
            return varmask
        
        
def getvarens(varf):
    for varens in ["ensemble", "ensembles", "M", "realization"]: #, "lev", "height"]:
        if varens in varf.variables.keys():
            return varens

def getdimens(varf, varname):
    for varens in ["ensemble", "ensembles", "M"]: #, "lev", "height"]:
        if varens in varf.variables[varname].dimensions:
            return varens     
        
def getdimlon(varf, varname):
    for varlon in ["x", "lon"]:
        if varlon in varf.variables[varname].dimensions:
            return varlon     
        
def getdimlat(varf, varname):
    for varlat in ["y", "lat"]:
        if varlat in varf.variables[varname].dimensions:
            return varlat    
        



def extract_array(varf, varname, nmon, lon_bnds, lat_bnds, level="all"):
    """
    extract a region from a netcdf file
    might have issue with certain grid (irregular or not splitted at greenwitch)
    varf: netcdf file open
    varname/ name of the variable to read
    lon_bnds: limits of longitude of the box to exctract
    lat_bnds: limits of latitude of the box to exctract
    level: the level to select within the file 
    """
    #print(level)
    varfvar=varf.variables[varname]
    varlat=getvarlat(varf)
    varlon=getvarlon(varf)
        
    lats = varf.variables[varlat][:] 
    lons = varf.variables[varlon][:]
    #print lons[0]
    #lat_inds = np.where((lats >= lat_bnds[0]) & (lats <= lat_bnds[1]))[0]
    #lon_inds = lon_index(lons, lon_bnds )
    
    lat1D, lon1D, lat_inds, lon_inds = lonlat_index(lats, lons, lat_bnds, lon_bnds)
    #print 
        #print(varf.variables[varname])
    try:
        unit=varf.variables[varname].units
    except:
        defaultdic={"tos":"K","ts":"K", "tauuo":"N m**-2", "tauu":"N m**-2"}
        unit=defaultdic.get(varname)
        print("warning: unit not found in the file, set to default: "+unit)
        
    #print(unit)
    offsetdic={"Celsius_scale":0, "K":-273.15, "mm/day":0, "m s-1":0, "m s**-1":0, "m/s":0, 
               "Kelvin_scale":-273.15, "N m**-2 s":0, "degC":0, "N m**-2":0, "N/m2":0, 
              "m s**-1":0, "m/s":0, "m":0, "Pa":0}
    scaledic={"Celsius_scale":1, "K":1, "mm/day":1,
              "Kelvin_scale":1, "N m**-2 s":1./21600, "degC":1, "N m**-2":1, "N/m2":1, 
              "m s**-1":86400*1000, "m/s":86400*1000,"m s-1":86400*1000, "m":1000, "Pa":1./100}
    if varname in ["uas", "vas"]:
        scaledic={"m s-1":1,"m s**-1":1, "m/s":1, "m s**-1":1}
   

    offset = offsetdic.get(unit)
    scale=scaledic.get(unit)
    
    if scale==None:
        scale=1
    if offset==None:
        offset=0
    #print(scale, offset)
    ndim=len(varf.variables[varname].shape)
    
    print(ndim)
    
    #print(varf.variables[varname].shape)
    #if str(varf.variables[varname].dimensions[1])==getdimlat(varf, varname):
    #    varfvar=np.swapaxes(vararray,1,2)
    #print varfvar.shape
    #print "lon inds", lon_inds
    #print len(lon_inds)
    if len(lon_inds)==2:
        if  ndim==5:
            if level!="all":
                vararrayW=varfvar[0:nmon,:,level,lat_inds[0]:(lat_inds[-1]+1),lon_inds[0][0]:(lon_inds[0][-1]+1)]
                vararrayE=varfvar[0:nmon,:,level,lat_inds[0]:(lat_inds[-1]+1),lon_inds[1][0]:(lon_inds[1][-1]+1)]
                #print unit,scaledic.get(unit),offsetdic.get(unit)
                vararray=np.concatenate((vararrayW, vararrayE), axis=3)*scale+offset
            else:
                vararrayW=varfvar[0:nmon,:,:,lat_inds[0]:(lat_inds[-1]+1),lon_inds[0][0]:(lon_inds[0][-1]+1)]
                vararrayE=varfvar[0:nmon,:,:,lat_inds[0]:(lat_inds[-1]+1),lon_inds[1][0]:(lon_inds[1][-1]+1)]
                #print unit,scaledic.get(unit),offsetdic.get(unit)
                vararray=np.concatenate((vararrayW, vararrayE), axis=4)*scale+offset
        elif ndim==4:
            if level!="all":
                vararrayW=varfvar[0:nmon,level,lat_inds[0]:(lat_inds[-1]+1),lon_inds[0][0]:(lon_inds[0][-1]+1)]
                vararrayE=varfvar[0:nmon,level,lat_inds[0]:(lat_inds[-1]+1),lon_inds[1][0]:(lon_inds[1][-1]+1)]
                #print unit,scaledic.get(unit),offsetdic.get(unit)
                vararray=np.concatenate((vararrayW, vararrayE), axis=2)*scale+offset
            else:
                vararrayW=varfvar[0:nmon,:,lat_inds[0]:(lat_inds[-1]+1),lon_inds[0][0]:(lon_inds[0][-1]+1)]
                vararrayE=varfvar[0:nmon,:,lat_inds[0]:(lat_inds[-1]+1),lon_inds[1][0]:(lon_inds[1][-1]+1)]
                #print unit,scaledic.get(unit),offsetdic.get(unit)
                vararray=np.concatenate((vararrayW, vararrayE), axis=3)*scale+offset
                
        elif ndim==3:
            vararrayW=varfvar[0:nmon,lat_inds[0]:(lat_inds[-1]+1),lon_inds[0][0]:(lon_inds[0][-1]+1)]
            vararrayE=varfvar[0:nmon,lat_inds[0]:(lat_inds[-1]+1),lon_inds[1][0]:(lon_inds[1][-1]+1)]
            vararray=np.concatenate((vararrayW, vararrayE), axis=2)*scale+offset
            
            
        if len(lons.shape)==1:    
            lons_regW=lons[lon_inds[0]]
            lons_regE=lons[lon_inds[1]]
            lons_reg=np.concatenate((lons_regW, lons_regE), axis=0)
            lats_reg=lats[lat_inds]
        elif len(lons.shape)==2:
            lonsW=lons[lat_inds[0]:(lat_inds[-1]+1),lon_inds[0][0]:(lon_inds[0][-1]+1)]
            lonsE=lons[lat_inds[0]:(lat_inds[-1]+1),lon_inds[1][0]:(lon_inds[1][-1]+1)]

            lons_reg=np.concatenate((lonsW, lonsE), axis=1)
            
            latsW=lats[lat_inds[0]:(lat_inds[-1]+1),lon_inds[0][0]:(lon_inds[0][-1]+1)]
            latsE=lats[lat_inds[0]:(lat_inds[-1]+1),lon_inds[1][0]:(lon_inds[1][-1]+1)]
            #print("I am here")
            #print(latsW.shape, latsE.shape)
            lats_reg=np.concatenate((latsW, latsE), axis=1)
            #print(lats_reg.shape)

    else:
                #print(lon_inds[0])
        if ndim==5:
            if level!="all":
                vararray=varfvar[0:nmon,:,level,lat_inds[0]:(lat_inds[-1]+1),
                        lon_inds[0][0]:(lon_inds[0][-1]+1)]*scale+offset
            else:
                vararray=varfvar[0:nmon,:,:,lat_inds[0]:(lat_inds[-1]+1),
                        lon_inds[0][0]:(lon_inds[0][-1]+1)]*scale+offset
        if ndim==4:
            if level!="all":
                vararray=varfvar[0:nmon,level,lat_inds[0]:(lat_inds[-1]+1),
                        lon_inds[0][0]:(lon_inds[0][-1]+1)]*scale+offset
            else:
                vararray=varfvar[0:nmon,:,lat_inds[0]:(lat_inds[-1]+1),
                        lon_inds[0][0]:(lon_inds[0][-1]+1)]*scale+offset
        elif ndim==3:
            vararray=varfvar[0:nmon,lat_inds[0]:(lat_inds[-1]+1),
                        lon_inds[0][0]:(lon_inds[0][-1]+1)]*scale+offset
        
        #print lons
        
        
        if len(lons.shape)==1:   
            print("I am here")
            print(len(lons.shape))
            lons_reg=lons[lon_inds[0]]
            lats_reg=lats[lat_inds]

        elif len(lons.shape)==2:
            lats_reg=lats[lat_inds[0]:(lat_inds[-1]+1),lon_inds[0][0]:(lon_inds[0][-1]+1)]
            lons_reg=lons[lat_inds[0]:(lat_inds[-1]+1),lon_inds[0][0]:(lon_inds[0][-1]+1)]


    #if lats[0]>lats[1]:
    #    lats=lats[::-1]
    #    if ndim==4:
    #        vararray=vararray[:,:,::-1,:]
    #    if ndim==3:
    #         vararray=vararray[:,::-1,:]
           
        #print("lat reversed")  
    print(vararray.shape)

    return(vararray, lats_reg, lons_reg) 



