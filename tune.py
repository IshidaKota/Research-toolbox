import numpy as np
from cmocean import cm
import matplotlib.pyplot as plt
from pyproj import Proj
import sys
import warnings
import pandas as pd
import math


class Fvcom2D:
    '''
    Postprocessing for FVCOM 4 netcdf output in its original coordinates
    '''
    def __init__(self, ncfile_path='tst_0001.nc', 
                 obcfile_path='tst_obc.dat', m_to_km=True, offset=False):
        '''
        Parameters
        ----------
        ncfile_path : str
            FVCOM output netcdf file
        obcfile_path : str
            FVCOM casename_obc.dat file
        m_to_km : bool
            = True if converting x and y axis units from m to km.
        offset : bool
           offset = True if the left-bottom corner is set to the origin.

        Returns
        -------
        Instance of class FVCOM2D
        '''

        self.fvcom_nc=ncfile_path
        if os.path.isfile(self.fvcom_nc):
            self.FVCOM = netCDF4.Dataset(self.fvcom_nc, 'r')
        else:
            print(f"ERROR: File {self.fvcom_nc} does not exit.")
            sys.exit()
        self.variables = self.FVCOM.variables
        print(f"Dictionaly keys of netCDF4.Dataset = {self.variables.keys()}")
        self.fobc=obcfile_path # = None if not exist
        if os.path.isfile(self.fobc):
            df = pd.read_csv(self.fobc, header=None, skiprows=1, delim_whitespace=True)
            ### -1 because index in FVCOM starts from 1 while from 0 in Python
            self.node_bc = df.iloc[:,1].values - 1
            print(f"Open boundary nodes = {self.node_bc}")
        else:
            print(f"ERROR: File {self.fobc} does not exit.")
            sys.exit()            
        self.x, self.y = self.variables['x'][:], self.variables['y'][:] ### at node
        self.xc, self.yc = self.variables['xc'][:], self.variables['yc'][:]   ### at cell center
        # Convert axis units from m to km
        if m_to_km:
            self.x /= 1000.0; self.y /= 1000.0; self.xc /= 1000.0; self.yc /= 1000.0
        # Set the left-bottom corner as the coordinate origin. 
        if offset:
            xmin = self.x.min(); ymin = self.y.min()
            self.x -= xmin; self.y -= ymin; self.xc -= xmin; self.yc -= ymin
        # Gets sigma coordinate values
        self.siglay, self.siglev = self.variables['siglay'][:], self.variables['siglev'][:]
        # Get time variables
        self.iint = self.variables['iint'][:]
        self.time = self.variables['time'][:]
        self.Itime = self.variables['Itime'][:]
        self.Itime2 = self.variables['Itime2'][:]
        # Gets bathymetry
        self.z = self.variables['h'][:]
        # Gets verts = [(x0,y0,h0),[x1,y1,h1],...,[xn,yn,hn]]
        verts = [(xi, yi, hi) for xi, yi, hi in zip(self.x.data, self.y.data, self.z.data)]
        self.verts = pd.DataFrame(verts, columns=['x', 'y', 'z'])
        # Gets connectivity array nv (3 node numbers of element i)
        nv = self.variables['nv'][:].T - 1
        self.triang = tri.Triangulation(self.x, self.y, triangles=nv)
        # Since node indexes in an element are defined clockwise in FVCOM,
        #     change them to counterclockwise, which is the matplotlib.tri specification.
        # FVCOMでは時計回りに定義されているので，matplotlib.triの仕様である反時計回りに変更する．
        # node番号自体は不変．
        self.nv=nv[:,::-1]
        self.tris = pd.DataFrame(self.nv, columns=['v0', 'v1', 'v2'])
        self.mesh = du.mesh(self.verts, self.tris) 
        self.trimesh = hv.TriMesh((self.tris, self.verts))
        # Element index of 3 neighbors of element i
        # 各セルには一般に3つのセルが隣接（境界セルはこの限りではない）
        self.nbe = np.array([[self.nv[n, j], self.nv[n, (j+2)%3]] \
            for n in range(len(self.triang.neighbors)) \
            for j in range(3) if self.triang.neighbors[n,j] == -1])

    @property
    def timestamp(self):
        '''
        Create time stamp list in ["D HH:MM:SS"]
        '''
        dayi = self.Itime
        hourf = self.Itime2 / 3600000
        houri = hourf.astype('int32')
        minf = (hourf - houri) * 60
        mini = minf.astype('int32')
        secf = (minf - mini) * 60
        seci = secf.astype('int32')
        return [f"{d} {h:02}:{m:02}:{s:02}" for d, h, m, s in zip(dayi, houri, mini, seci)]
    
#ここから    
#get variables
f = '../examples/Tokyo2_0001.nc'
obcfile = 'C:/Users/ishid/fvcominputs/input/TokyoBay_obc.dat'
fvcom = Fvcom2D(ncfile_path=ncfile, obcfile_path=obcfile, m_to_km=True, offset=False)
sal = nc.variables['salinity']
temp = nc.variables['temprature']

#get station node
stn ={'chiba1buoy':(139.9542517,35.53703833),'chibaharo':(140.0233033,35.61095833),'kawasaki':(139.8340267,35.49019),'urayasu':(139.9417417,35.640085)}

def find_closest_node(stn):
    index = {}
    p = Proj(proj='utm',zone=54,ellps='WGS84', preserve_units=False)
    
    for key in stn.keys():
        st_m = p(stn[key][0],stn[key][1]) #lon,lat→utm(m)
        st_km= st_m[0]/1000,st_m[1]/1000 #utm(m)→utm(km)
        #それぞれの距離を計算、listで保存
        dist = [math.sqrt((fvcom.x[i]-st_km[0])**2 +(fvcom.y[i]-st_km[1])**2) for i in range(len(fvcom.x))]
        #listの中から最小値を見つける
        res = dist.index(min(dist))
        min_dist = min(dist)
        #辞書に追加
        index[key] = res
        print(f"Node {res} is {min_dist}km to the station:{key}")
        
    return index

stn_node = find_closest_node(stn)
print(stn_node)

def calc_cc(x,y):
    s1 = pd.Series(x)
    s2 = pd.Series(y)
   # print(s1,s2)
    res = s2.corr(s1)
    return res

#get observation data to plot each staion temprature and plot
def get_temp_cc(stn_node):
    cc = {}
    rmse = {}
    for key in stn_node:
        index = stn_node[key]
        #read Mpos_s
        f = 'C:/Users/ishid/Github/Puppeteer/src/csv/Mpos_S_' + key +'_2019.csv'
        df = pd.read_csv(f)
        
        df_tp_1 = df[df['lev']==1]#明日までに水深を考慮して計算する。sigma layerと深さの実数値を変換して、観測データのほうを内挿補間して比較する
        #df_tp_1 = df_tp_1.interpolate()
        df_tp_1 = df_tp_1.reset_index(drop=True)
        df_tp_1.to_csv('C:/Users/ishid/Github/Puppeteer/src/csv/Mpos_S_' + key +'interp_2019.csv')
        
        #plot
        fig = plt.figure(figsize=(50,4))
        ax = fig.add_subplot(1,1,1,xlabel = 'hours since 2019/1/1 0:00:00(UTC)',ylabel='temprature degree')
       
        ax.plot(temp[:,0, index],linewidth=2,label='simulation') #float64 temp(time, siglay, node)
        ax.plot(df_tp_1['tp'],label='observation')
        
        fig.suptitle('temprature@{},min mesh size800m'.format(key))
        ax.legend()
        plt.grid()
        plt.plot()
        #plt.show()
        fig.savefig('temp_'+key+'.png')
        
        #calc RMSE
        n = 24 #とりあえずバグを起こさないように(観測値とモデルで長さがあっていない)
        res = 0
        for i in range(n):
            res += (temp[i,0,index] -df_tp_1.iloc[i,6])**2 #temp:6
            res = res/n
            res = math.sqrt(res)
            
        rmse[key] = res
        
        #calc cc
        res = calc_cc(temp[:,0, index],df_tp_1['tp'])
        cc[key] = res
    return rmse,cc

cc = get_temp_cc(stn_node)
print(cc)
    
    
    
    
    
    
    
    
    

    
