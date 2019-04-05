import sys
import os
# import Nio
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import numba
from numba import jit
import dask
from dask import delayed
import psutil
from rasterio.plot import show

def main():
    inpfname = "/home/gabriel/transicao/doutorado/rochedo/OTIMIZAGRO_Marcos_Costa/WEG/weg_comp_2012-2050.nc"
    inpvname = "landuse"
    codfname = "/home/gabriel/transicao/doutorado/rochedo/OTIMIZAGRO_Marcos_Costa/biomas.nc"
    codvname = "Band1"

    cods = [1,2,3,4,5,6]
    clas = (list(range(0,40)))

    maxcod = 6
    maxcla = 40


    inpvarname = "PRECT"
    outvarname = "prec"

    finp = xr.open_dataset(inpfname, decode_times=False)
    inparr = finp[inpvname]
    times = np.round(inparr['time'])

    fcod = xr.open_dataset(codfname)
    codarr = fcod[codvname]

    year = times[0]

    yinparr = inparr.sel(time = year, method = 'nearest')

    yinparr = yinparr.sel(lat = slice(-10.0,0.0), lon = slice(-70.,-60.))
    codarr = codarr.sel(lat = slice(-10.0,0.0), lon = slice(-70.,-60.))
    # yinparr = yinparr.sel(lat = slice(-2.0,0.0), lon = slice(-62.,-60.))
    # codarr = codarr.sel(lat = slice(-2.0,0.0), lon = slice(-62.,-60.))
    # show(yinparr)
    # show(codarr)

    cod = 5

    # vec = yinparr.values[10,:]
    # # print(sum1d(vec)
    # mat = yinparr.values[10:20,10:20]
    # codmat = codarr.values[10:20,10:20]
    # print(mat)
    # print(yinparr.dtype)
    # print(sum2d(mat))
    # print(sum2d(yinparr.values))

    # outmat = np.zeros((len(cods),len(clas)), dtype = inparr.dtype)
    outmat = np.zeros((maxcod+1,maxcla+1), dtype = np.float64)

    chunksize = 100
    # print(countclasses(yinparr.values.astype(int),outmat))
    array = dask.array.from_array(yinparr.values.astype(int), chunks = chunksize)
    codarray = dask.array.from_array(codarr.values.astype(int), chunks = chunksize)
    # array = yinparr.values.astype(int)
    # codarray = codarr.values.astype(int)
    print(sys.getsizeof(array))
    print(sys.getsizeof(codarray))
    # exit()
    # filloutmat = countclasses_cod(array, codarray, outmat)
    # filloutmat = countclassesc_cod(mat.astype(int),codmat.astype(int), outmat)
    # print(np.isfinite(codarr.values[0,0]))
    # lala = filloutmat.compute()

    # icods = {cods[i]:i for i in range(len(cods))}
    # iclas = {clas[i]:i for i in range(len(clas))}
    # print(iclas)

@delayed
# @jit(nopython = True, parallel = True)
@jit(nopython = True)
def countclasses_cod(array,codarray,outmat):
    # for i in numba.prange(array.shape[0]):
        # for j in numba.prange(array.shape[1]):
    for i in range(array.shape[0]):
        for j in range(array.shape[1]):
            if np.isfinite(codarray[i,j]) & np.isfinite(array[i,j]):
                # print(outmat[codarray[i,j],array[i,j]])
                outmat[codarray[i,j],array[i,j]] += 1
    # print(np.max(outmat))
    return outmat

@jit(nopython = True)
def countclasses(array,outmat):
    for i in range(array.shape[0]):
        for j in range(array.shape[1]):
            # sum += array[i,j]
            outmat[0,array[i,j]] += 1.0
    return outmat


@jit(nopython = True)
def countclasses(array,outmat):
    for i in range(array.shape[0]):
        for j in range(array.shape[1]):
            # sum += array[i,j]
            outmat[0,array[i,j]] += 1.0
    return outmat


@jit(nopython = True)
def sum2d(array):
    sum = 0.0
    for i in range(array.shape[0]):
        for j in range(array.shape[1]):
            sum += array[i,j]
    return sum


@jit
def sum1d(array):
    sum = 0.0
    for i in range(array.shape[0]):
        sum += array[i]
    return sum

    # print(yinparr.lat)
    # print(codarr.lat)
    #
    # mskyinparr = yinparr
    # mskyinparr = mskyinparr.where(codarr == cod)
    # show(mskyinparr)
    #
    # print(np.sum(mskyinparr))
    # # cod = 5
    #
    # # MemoryError
    # indcodes = {cod:np.where(codarr == cod) for cod in cods}
    #
    # print(np.sum(yinparr.values[indcodes[cod]]))



    # print(indcodes)

    # mskinparr = inparr.
    #
    # inparr.sel(time = year, method = 'nearest').where(codarr == cod)

    #inparr.sel(time = 2012, method = 'nearest')

    # finp.to_netcdf(outfname, mode = "w", format = "NETCDF3_64BIT")
    # finp[inpvarname].to_netcdf(outfname, mode = "w", format = "NETCDF3_64BIT")




if __name__ == "__main__":
    main()
