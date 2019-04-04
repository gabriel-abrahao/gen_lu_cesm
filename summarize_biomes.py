import os
# import Nio
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import numba
from numba import jit
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
    # show(yinparr)

    cod = 5

    # vec = yinparr.values[10,:]
    # print(sum1d(vec))
    mat = yinparr.values[10:20,10:20]
    print(mat)
    print(yinparr.dtype)
    print(sum2d(mat))
    # print(sum2d(yinparr.values))

    # outmat = np.zeros((len(cods),len(clas)), dtype = inparr.dtype)
    outmat = np.zeros((maxcod+1,maxcla+1), dtype = inparr.dtype)
    print(countclasses(yinparr.values.astype(int),outmat))

    # icods = {cods[i]:i for i in range(len(cods))}
    # iclas = {clas[i]:i for i in range(len(clas))}
    # print(iclas)

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
