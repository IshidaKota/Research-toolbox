## Load an FVCOM model output and plot surface
from PyFVCOM.read import FileReader
from PyFVCOM.plot import Time
from PyFVCOM.grid import unstructured_grid_depths
import cmocean as cm 
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from pyproj import Proj
import sys



def find_closest_node(stn):
    """
    moniteringpostから最近傍の場所を抽出
    """
    index = {}
    p = Proj(proj='utm',zone=54,ellps='WGS84', preserve_units=False)
    
    for key in stn.keys():
        st_m = p(stn[key][0],stn[key][1]) #lon,lat→utm(m)
        st_km= st_m[0]/1000,st_m[1]/1000 #utm(m)→utm(km)
        dist = [math.sqrt((fvcom.x[i]-st_km[0])**2 +(fvcom.y[i]-st_km[1])**2) for i in range(len(fvcom.x))] #全nodeからの距離を計算、listで保存
        #listの中から最小値を見つける
        res ,min_dist= dist.index(min(dist)),min(dist) #indexを取得
        index[key] = res
        print(f"Node {res} is {min_dist}km to the station:{key}")
    return index

#stn ={'chiba1buoy':(139.9542517,35.53703833),'chibaharo':(140.0233033,35.61095833),'kawasaki':(139.8340267,35.49019),'urayasu':(139.9417417,35.640085)}
#stn_node = find_closest_node(stn)
stn = {'chiba1buoy':136,'chibaharo':46,'kawasaki':222,'urayasu':51}

#fnames= [sys.argv[i] for i in range(1,len(sys.argv)-1)]##only before _0001.nc
fnames = [sys.argv[1]]
for fname in fnames:
    ncfile = "../run/"+fname + "_0001.nc"
    print(ncfile)
    for place in stn:
        index = stn[place]
    #gauge = (140.05, 35.57)  # a sample (lon, lat) position for Estuary example
    #index = fvcom.closest_node(gauge).item()  ### ndarray -> scalar
    # print(index) # 295
        fvcom = FileReader(ncfile, variables=['zeta', 'temp'], dims={'time': slice(0,8700), 'node': index})
        time = Time(fvcom, figsize=(35,12), cmap = 'jet',cb_label='{} ({})'.format(fvcom.atts.temp.long_name,
                                                                    fvcom.atts.temp.units))
        z = unstructured_grid_depths(fvcom.grid.h, fvcom.data.zeta, fvcom.grid.siglay)
        ## print("z.shape", z.shape)
        ## fill_seabed makes the part of the plot below the seabed grey.
        ## We need to squeeze the data array since we've only extracted a single position.
        ## print(fvcom.grid.h.shape)
        ## print(fvcom.grid.h)
        time.plot_surface(z, fvcom.data.temp, fill_seabed=True, h_offset=0, h_min=fvcom.grid.h,norm=Normalize(vmin=5, vmax=30))
        time.axes.set_ylabel('Depth (m)')
        time.figure.savefig('./png/contour/'+fname+'_'+place+'.png')
        print('savefig,path=./png/contour/'+fname+'_'+place+'.png')