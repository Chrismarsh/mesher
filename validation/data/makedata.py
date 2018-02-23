import rasterio
import os
import numpy as np
import matplotlib.pyplot as plt

def flat(input_dem, output_dem, elevation):
    # Open dem_flat to read in info
    dem_flat_in = rasterio.open(input_dem)
    dem     = dem_flat_in.read(1)

    ###### Create flat dem
    dem_flat = np.array(dem, copy=True, dtype=np.float64)
    dem_flat.fill(elevation)

    # Write out to geotiff
    with rasterio.open(output_dem, 'w', driver='GTiff', height=dem_flat_in.shape[0],
                       width=dem_flat_in.shape[1], count=1, dtype=rasterio.float64,
                       crs=dem_flat_in.crs, transform=dem_flat_in.transform) as dst: #'+proj=latlong'
        dst.write(dem_flat, 1)

    return rasterio.open(output_dem)

def vegpatch_flat(dem_flat,output_dem):
    dem_veg = np.array(dem_flat.read(1), copy=True, dtype=np.float64)
    x = int(dem_veg.shape[0]/2.)+50
    y = int(dem_veg.shape[1]/2.)+50

    dem_veg[:] = 1
    ps = 500 # 1/2 patch size
    dem_veg[x-ps:x+ps,y-ps:y+ps] = 100.0

    with rasterio.open(output_dem, 'w', driver='GTiff', height=dem_flat.shape[0],
                       width=dem_flat.shape[1], count=1, dtype=rasterio.float64,
                       crs=dem_flat.crs, transform=dem_flat.transform) as dst: #'+proj=latlong'
        dst.write(dem_veg, 1)

    return rasterio.open(output_dem)

def gausshill(dem_flat,output_dem):

    dem_gauss = np.array(dem_flat.read(1), copy=True, dtype=np.float64)
    hm = 1000 # Mountain height (m)

    def makeGaussian(size, fwhm=3, center=None):
        """ Make a square gaussian kernel.
        size is the length of a side of the square
        fwhm is full-width-half-maximum, which
        can be thought of as an effective radius.
        """
        x = np.arange(0, size, 1, float)
        y = x[:, np.newaxis]

        if center is None:
            x0 = y0 = size // 2
        else:
            x0 = center[0]
            y0 = center[1]

        return np.exp(-4 * np.log(2) * ((x - x0) ** 2 + (y - y0) ** 2) / fwhm ** 2)

    # Make Gaussian hill (height 0-1)
    dem_gauss = makeGaussian(np.max(dem_flat.shape), fwhm=np.max(dem_flat.shape)/4., center=(dem_flat.shape[0]/2,dem_flat.shape[1]/2))

    # Scale to real elevation range (500m for example)
    dem_gauss = dem_gauss * hm

    # Extract only cells from orig dem shape
    dem_gauss_out = dem_gauss#[:377, :518]

    # Add in base elevation
    dem_gauss_out += (dem_gauss - hm)

    # Write out to geotiff
    print(dem_gauss_out.shape)
    with rasterio.open(output_dem, 'w', driver='GTiff', height=dem_gauss_out.shape[0],
                       width=dem_gauss_out.shape[1], count=1, dtype=rasterio.float64,
                       crs=dem_flat.crs, transform=dem_flat.transform) as dst: #'+proj=latlong'
        dst.write(dem_gauss_out, 1)

    return rasterio.open(output_dem)


def idealridge(dem_flat,output_dem):

    dem_ridge = np.array(dem_flat.read(1), copy=True, dtype=np.float64)
    dxdy = flat_dem.get_transform()[1] # resolution of dem_flat dem, assume square cells
    hm = 300 # Mountain height (m)
    # xa = 1200/dxdy # Half width (m)
    xa = 600/dxdy # Half width (m)
    icm = dem_ridge.shape[0] / 2

    # Create idealized ridge
    for (x,y), value in np.ndenumerate(dem_ridge):
        dem_ridge[x,y] = hm/(1.+(float(x-icm)/xa)**2)

    # Add in base elevation
    dem_ridge += (dem_ridge - hm)

    # Write out to geotiff
    print(dem_ridge.shape)
    with rasterio.open(output_dem_ridge, 'w', driver='GTiff', height=dem_ridge.shape[0],
                       width=dem_ridge.shape[1], count=1, dtype=rasterio.float64,
                       crs=dem_flat.crs, transform=dem_flat.transform) as dst: #'+proj=latlong'
        dst.write(dem_ridge, 1)

    return rasterio.open(output_dem)


main_dir = ''
input_dem  = os.path.join(main_dir,'granger1m.tif') # Uses this DEMS resolution and projection info
output_dem_flat = os.path.join(main_dir, 'ideal_flat.tif')
output_dem_ridge = os.path.join(main_dir, 'ideal_ridge.tif')
output_dem_guass = os.path.join(main_dir, 'ideal_gaussianHill.tif')
output_dem_veg = os.path.join(main_dir, 'ideal_flatveg.tif')

flat_dem = flat(input_dem,output_dem_flat,1000.)
# plt.figure()
# plt.imshow(flat_dem.read(1))
# plt.colorbar()

ridge = idealridge(flat_dem,output_dem_ridge)
plt.figure()
plt.imshow(ridge.read(1))
plt.colorbar()

gauss = gausshill(flat_dem,output_dem_guass)
plt.figure()
plt.imshow(gauss.read(1))
plt.colorbar()

veg = vegpatch_flat(flat_dem,output_dem_veg)
plt.figure()
plt.imshow(veg.read(1))
plt.colorbar()

# plt.show()