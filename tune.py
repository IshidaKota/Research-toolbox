
"""
Usage:
python tune.py input.file 

without extention
"""

import pandas as pd
import numpy as np
import math
import sys
import os
from functools import partial, partialmethod
import pyproj
from pyproj import Proj
### import xarray as xr
import netCDF4
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import holoviews as hv
from holoviews import opts
import datashader.utils as du
from holoviews.operation.datashader import datashade, shade, dynspread, rasterize,\
    spread, aggregate, regrid
from holoviews.operation import decimate
import holoviews.plotting.mpl
import geoviews as gv
from bokeh import models
import subprocess
from bokeh.io import export_png
from scipy import interpolate
#hv.extension('bokeh', 'matplotlib')


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
        
    def time_list(self):
        Time_list = []
        Times = self.variables['Times']
        for i in range(len(Times)):
            time =Times[i].data
            bytes_time = b''
            for j in range(19):
                bytes_time = bytes_time + time[j]
            time2 = bytes_time.decode('utf-8')
            time2 = time2.replace('T',' ')
            Time_list.append(time2)
        return Time_list

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

    def plot_mesh(self, x_range=None, y_range=None, **kwargs):
        '''
        Plot mesh
        
        Parameters
        ----------
        x_range : tuple
            Plot range of x axis (xmin, xmax)
        y_range : tuple
            Plot range of y axis (ymin, ymax)
        '''

        # Gets **kwargs with default values
        #width=kwargs.get('width', 300)
        #height=kwargs.get('height', 300)
        line_color=kwargs.get('line_color', 'blue')
        line_width=kwargs.get('line_width', 0.1)

        # jsasaki debug
        #p = self.trimesh.edgepaths.options(width=width, height=height, 
        #                                   line_width=line_width, line_color=line_color)
        p = self.trimesh.edgepaths.options(line_width=line_width, line_color=line_color)

        # Set plot range
        if x_range is None and y_range is None:
            return p
        elif y_range is None:
            return p.redim.range(x=x_range)
        elif x_range is None:
            return p.redim.range(y=y_range)
        else:
            return p.redim.range(x=x_range, y=y_range)
    
    def plot_coastline(self, x_range=None, y_range=None, **kwargs):
        '''
        Plot coastline. May take much time for interactive plotting.

        Parameters
        ----------
        x_range : tuple
            Plot range in x-axis (xmin, xmax)
        y_range : tuple
            Plot range in y-aixs (ymin, ymax)
        '''

        # Gets **kwargs with default values
        #width=kwargs.get('width', 300)
        #height=kwargs.get('height', 300)
        color=kwargs.get('color', 'red')
        line_width=kwargs.get('line_width', 0.5)

        # paths = [[(x0s, y0s), (x0e, y0e)], ..., [(xns, yns), (xne, yne)]]
        paths = [[(self.x[self.nbe[m,:]][0], self.y[self.nbe[m,:]][0]), \
                  (self.x[self.nbe[m,:]][1], self.y[self.nbe[m,:]][1])] \
                 for m in range(len(self.nbe))]

        # jsasaki debug
        #p = hv.Path(paths).options(width=width, height=height, 
        #                           color=color, line_width=line_width)
        p = hv.Path(paths).options(color=color, line_width=line_width)

        if x_range is None and y_range is None:
            return p
        elif y_range is None:
            return p.redim.range(x=x_range)
        elif x_range is None:
            return p.redim.range(y=y_range)
        else:
            return p.redim.range(x=x_range, y=y_range)

    def get_val(self, item, time=None, sigma=None):
        '''
        Gets item values dependent considering its shape with specifying
        output time step (if time sereis data) and sigma level.

        Parameters
        ----------
        item : str
            Name of item in dictionary keys
        time : int
            Time series index
        sigma : int
            Layer number in sigma coordinates

        Returns
        -------
        pandas.DataFrame
        '''

        # Check whether item exists among the keys
        if not item in self.variables.keys():
            print(f"Error: Item {item} does not exit in keys.")
            sys.exit()
        val = self.FVCOM.variables[item]

        # Check the item dimensions
        if len(val.shape) == 1 : # 1-D item
            scalar = val
        elif len(val.shape) >= 2: # 2-D or 3-D item
            if time < 0 or time > val.shape[0]:
                print(f"ERROR: Time index = {time} is out of range.")
                print(f"Time index should be between 0 and {val.shape[0]-1}.")
                sys.exit()
            if len(val.shape) == 2: ### 2-D item
                scalar = val[time]
            elif len(val.shape) == 3: ### 3-D item
                if sigma < 0 or sigma > val.shape[1]:
                    print(f"ERROR: Sigma layer number = {sigma} is out of range.")
                    print(f"Sigma layer number should be between 0 and {val.shape[1]-1}.")
                    sys.exit()
                else:
                    scalar = val[time][sigma]
            else:
                print("ERROR: the shape is incorrect.")
                sys.exit()

        verts = [(xi, yi, hi) for xi, yi, hi in zip(self.x.data, self.y.data, scalar)]
        return pd.DataFrame(verts, columns=['x', 'y', item])

    def __plot_2dh(self, item, time, sigma):
        '''
        Internal method for converting pd.DataFrame of item to hv.TriMesh

        Parameters
        ----------
        item : str
        time : int
            Time index to be read
        sigma :int
            Sigma level index to be read
        Returns
        -------
        HoloViews.TriMesh
        '''

        # Reading pd.DataFrame of item using get.val() method 
        if len(self.variables[item].shape) == 1:   ### 1-D item (e.g.: depth) 
            df_item = self.get_val(item)
        elif len(self.variables[item].shape) == 2: ### 2-D item (e.g.: zeta)
            df_item = self.get_val(item, time=time)
        elif len(self.variables[item].shape) == 3: ### 3-D item (e.g.: salainity)
            df_item = self.get_val(item, time=time, sigma=sigma)
        else:
            print('ERROR: The shape of item is incorrect in self.__plot_2dh')
            sys.exit()
        return hv.TriMesh((self.tris, hv.Points(df_item)))

    def dmap(self, item, x_range=None, y_range=None, z_range=None, timestamp_location=None):
        '''
        Interactive plotting using DynamicMap
        Internal functions need to be deifned in order to be an argument of 
            functools.partial().
        Tips: Prepare an internal function,
            and set time and sigma at the end of the arguments in partial.
        内部関数を用意すること，partialの引数の最後をtimeとsigmaにするのがポイント
        
        Parameters
        ----------
        self : FVCOM2D instance
        item : str
        x_range : tuple
        y_range : tuple, optional
        z_range : tuple, optional
        timestamp_location : tuple, optional
        '''

        def item3d_plot(item, timestamp_location, sigma, time):
            '''
            Internal function for plotting 3D item using functools.partial().
            partialの引数とするため，内部関数とする.

            Parameters
            ----------
            item : str
            timestamp_location : tuple
            sigma : int
                sigma level
            time : int
                output timestep of time series data

            Returns
            -------
            hv.TriMesh * self.plot_timestamp
            '''

            item = self.get_val(item, time=time, sigma=sigma)
            timestamp_txt = self.plot_timestamp(time, timestamp_location)
            return hv.TriMesh((self.tris, hv.Points(item))) * timestamp_txt

        def item2d_plot(item, timestamp_location, time):
            '''
            Internal function for plotting 2D item (e.g. zeta) using functools.partial().
            
            Parameters
            ----------
            item : str
            timestamp_location : tuple
            time : int
            
            Returns
            -------
            hv.TriMesh * self.plot_timestamp
            '''

            item = self.get_val(item, time=time)
            timestamp_txt = self.plot_timestamp(time, timestamp_location)
            return hv.TriMesh((self.tris, hv.Points(item))) * timestamp_txt        
        
        if len(self.variables[item].shape) == 1:   ### 1-D item (ex: depth) 
            print(f"ERROR: 1-D item of {item} is not supported in FVCOM2D.dmap")
            sys.exit()
        elif len(self.variables[item].shape) == 2: ### 2-D item (ex: zeta)
            p = partial(item2d_plot, item, timestamp_location)
            dmap = hv.DynamicMap(p, kdims=['time']).redim.values(time=range(len(self.time)))
        elif len(self.variables[item].shape) == 3: ### 3-D item (ex: salainity)
            p = partial(item3d_plot, item, timestamp_location)
            dmap = hv.DynamicMap(p, kdims=['sigma', 'time']).redim.values( \
                sigma=range(len(self.siglay))).redim.values(time=range(len(self.time)))
        else:
            print('ERROR: The shape of the item is incorrect in FVCOM2D.dmap')
            sys.exit()

        if x_range is None:
            x_range = (self.x.min(), self.x.max())
        if y_range is None:
            y_range = (self.y.min(), self.y.max())
        if z_range is None:
            return rasterize(dmap, x_range=x_range, y_range=y_range)
        else:
            return rasterize(dmap.redim.range(**{item : z_range}), \
                             x_range=x_range, y_range=y_range)

    def splot_2dh(self, item, time=0, sigma=0, x_range=None, y_range=None,
                  z_range=None, timestamp_location=None, **kwargs):
        '''
        Static horizontal 2D plot (one sigma and time frame)
        静的な平面2次元プロット．InteractiveにはDynamicMapを用いたdmap()を使用する
        '''
        if x_range is None:
            x_range = (self.x.min(), self.x.max())
        if y_range is None:
            y_range = (self.y.min(), self.y.max())

        p = self.__plot_2dh(item=item, time=time, sigma=sigma)
        if z_range is None or z_range == False:
            p = rasterize(p, x_range=x_range, y_range=y_range)
        else:
            p = rasterize(p.redim.range(**{item : z_range}), x_range=x_range, y_range=y_range)

        timestamp_txt = self.plot_timestamp(time, timestamp_location)
        return p * timestamp_txt

    def plot_point(self, coords, color='red', marker='+', size=10, **kwargs):
        p = hv.Points(coords, **kwargs)
        return p.opts(color=color, marker=marker, size=size)
        
    def plot_timestamp(self, time=0, timestamp_location=None, **kwargs):
        '''
        timestamp_location: tuple (x, y)
        Currently setting the location using graphic space unsupported in Holoviews
        スクリーン座標系で位置が設定できるとよいが，現在はHoloviewsの仕様で不可
        timestamp_location=Noneのときは空文字を出力することで汎用化した
        '''
        if timestamp_location is None:
            return hv.Text(0, 0, "", **kwargs)
        else:
            timestamp = self.timestamp[time]
            return hv.Text(timestamp_location[0], timestamp_location[1], timestamp, **kwargs)

    def make_frame(self, item, title=None, time=None, sigma=None,
                   plot_mesh=True, plot_coastline=True,
                   x_range=None, y_range=None, z_range=None,
                   timestamp_location=None):
        '''
        Make a frame at one time step [and one sigma layer (only for 3D)]
        1つのframe出力．アニメーション作成に有効
        '''
        ptri = self.splot_2dh(item=item,time=time, sigma=sigma, colorbar_title=title,
                              x_range=x_range, y_range=y_range,
                              z_range=z_range, timestamp_location=timestamp_location)
        pmesh = self.plot_mesh(x_range=x_range, y_range=y_range)
        pcoastline = self.plot_coastline(x_range=x_range, y_range=y_range)
        #timestamp_txt = self.plot_timestamp(time, timestamp_location)

        if plot_mesh and plot_coastline:
            return ptri * pmesh * pcoastline
        elif plot_mesh:
            return ptri * pmesh
        elif plot_coastline:
            return ptri * pcoastline
        else:
            return ptri
        
        
        
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
        index[key] = res + 1 #+1 :nodeは1始まり、listは0始まり
        print(f"Node {res} is \
         {min_dist}km to the station:{key},station is at {stn[key]}")
        
    return index

def calc_cc(x,y):
    s1 = pd.Series(x)
    s2 = pd.Series(y)
   # print(s1,s2)
    res = s2.corr(s1)
    return res
    
def interp(df,Depth,var_name):
    df = df.dropna(subset=['tp','sl'])      #nanを落とす
    #df[8173*24:8174*24]
    dates = df['datetime'].unique()
    levmax = max(df['lev'])
    reg = {}
    for date in dates:
        mask = (df['datetime'] == date)
        tmp = df[mask]
        x  = np.array(tmp['depth'])
        if var_name == 'temperature':
            y  = np.array(tmp['tp'])
        else:
            y  = np.array(tmp['sl'])
            

        if x.size >= levmax/2:
            f = interpolate.interp1d(x, y,kind='linear',axis=-1,copy=True,bounds_error=None, \
                                       fill_value='extrapolate',assume_sorted=False)
            tp_new = f(np.array(Depth))
            reg[date] = tp_new
        else:
            reg[date] = np.array([0 for i in range(len(Depth))])
    return reg    



def plotter(stn_node,var_name): 
    if var_name == 'temperature':
        var = fvcom.variables['temp']
    elif var_name == 'salinity':
        var = fvcom.variables['salinity']
    else:
        raise ValueError(f"invalid variable '{var_name}'")

    cc = {}
    rmse = {}
    Time_list = fvcom.time_list()
    for key in stn_node:
        index = stn_node[key]
        #read Mpos_S
        #f = '/work/gy29/y29007/Github/data/csv/Mpos_S_' + key +'_2020.csv'
        #df = pd.read_csv(f,usecols=['depth','tp','sl','datetime','lev'])
        
        regress = regress_dict_tp[key] #dictのdict形式でまとまっている
        sigma =[2,10,18]
        for k in sigma:
            res = np.zeros(len(date_ary))
            i = 0
            for date in date_ary.values:
                if date in regress.keys(): #空データでなければ
                    if regress[date].all() == 0:
                        res[i] = np.nan
                    else:
                        res[i] = regress[date][k]
                if res[i] == 0:
                        res[i] = np.nan
                #debug_temperature
                if var_name == 'temperature':
                    if abs(res[i])>40:
                        print(date,res[i],i)
                        res[i] = np.nan
                #debug_salinity
                else:
                    if res[i]>50 or res[i]<0:
                        print(date,res[i],i)
                        res[i] = np.nan
                i+=1
            
    #plot
            fig = plt.figure(figsize=(10,4))
            ax = fig.add_subplot(1,1,1,xlabel = 'days since 2020/1/1 0:00:00(UTC)',ylabel='{}'.format(var_name))
            ax.plot(var[:,k, index],linewidth=2,label='simulation') #float64 temp(time, siglay, node)
            ax.plot(res[:],linewidth=2,label='observation')
            fig.suptitle('{}@{},sigma = {}/20'.format(var_name,key,k))
            ax.legend()
            plt.grid()
            plt.ylim(0,35)
            plt.xticks([24*30*i for i in range(12)], [30*i for i in range(12)])
            plt.plot()
            os.makedirs('./png/'+sys.argv[1],exist_ok=True)
            fig.savefig('./png/'+sys.argv[1]+'/'+var_name+'_'+key+'_'+str(k)+'_'+sys.argv[1]+'.png')
            print('save figure@sigma={},stn={}'.format(k,key))
            
    #calc RMSE and CC
            #RMSE
            X = []
            Y = []
            diff = 0
            for i in range(len(Time_list)):
                if Time_list[i] in regress.keys() and res[i] != np.nan:
                    X.append(var[i,k,index].item())
                    Y.append(res[i])
                    diff += (var[i,k,index].item()-res[i])**2
                    if res[i] == np.nan:
                        print(i)
            print(f"Var is {diff}") #debug
            Rmse = math.sqrt(diff/len(Y)) #root 1/n

            #CC
            Cc = calc_cc(X,Y)

            cc[key] = Cc
            rmse[key] = Rmse
    return rmse,cc






f = '../run/'+sys.argv[1]+ '_0001.nc'
####################################USER DIFINITION#########################################################################
date_index = pd.date_range("2020-01-01 00:00:00", periods=24*366, freq="H") #specify start date
date_ary = date_index.to_series().dt.strftime("%Y-%m-%d %H:%M:%S")
obcfile = '../input_testcase2/TokyoBay_obc.dat' #specify where obcfile is
####################################END OF USER DIFINITION##################################################################

fvcom = Fvcom2D(ncfile_path=f, obcfile_path=obcfile, m_to_km=True, offset=False)
temp = fvcom.variables['temp']
sal = fvcom.variables['salinity']
stn ={'chiba1buoy':(139.9542517,35.53703833),'chibaharo':(140.0233033,35.61095833),'kawasaki':(139.8340267,35.49019),'urayasu':(139.9417417,35.640085)}
stn_node = find_closest_node(stn)

regress_dict_tp = {}
for key in stn_node:
    index = stn_node[key]
    #f = '/work/gy29/y29007/Github/data/csv/Mpos_S_' + key +'_2020.csv' #OBCX
    f = './data/Mpos_S_' + key +'_2020.csv'
    df = pd.read_csv(f)
    
    siglay = -fvcom.siglay[:,index].data
    Depth = [fvcom.variables['h'][index]*siglay[i] for i in range(len(siglay))]

    regress = interp(df,Depth,'temperature')
    regress_dict_tp[key] = regress

    
cc = plotter(stn_node,'temperature')
print(cc)
with open('result_'+sys.argv[1]+'.txt','w') as f:
    print(f"{sys.argv[1]} ,cc and RMSE are {cc}",file=f)
    print(find_closest_node(stn),file=f)
    
    
    
    
    
    
    
    
    
"""   
##########################################MAKE 3D PLOT HTML##################################################################
    
    
    
# Select item to be plotted.
item='temp'; title = "temprature"; z_range=(15,35)

# Set x, y, and z ranges. = None for auto setting (you may change after plotting).
x_range=None
y_range=None

# Set timestamp_location. = None to deactivate
timestamp_location=(270, 148)

# Create color map.
dmap_ds=fvcom.dmap(item=item, x_range=x_range, y_range=y_range, z_range=z_range, \
                   timestamp_location=timestamp_location)
# Create mesh.
mesh = fvcom.plot_mesh(line_color='black', line_width=0.1)

# Create coastline.
coastline = fvcom.plot_coastline(color='black', line_width=1)

# Overlay as you like
#p = dmap_ds * mesh * coastline

# Apply customization options.
# Customize colorbar.
colorbar_opts=dict(title=title, title_text_font_size="16pt", 
                   title_text_baseline = "bottom",
                   title_standoff =10, title_text_font_style = "normal",
                   major_label_text_font_size = "14pt", major_label_text_align="left",
                   label_standoff=4, width=15)
# Customize Image plot.
p_opts = opts.Image(tools=['hover'], cmap='rainbow', colorbar_position='right',
                    frame_width=400, aspect='equal', colorbar=True, colorbar_opts=colorbar_opts)
#p.redim.label(x='X (km)', y='Y (km)').opts(p_opts)
dmap_ds = dmap_ds.redim.label(x='X (km)', y='Y (km)').opts(p_opts)
# Overlay as you like
p = dmap_ds * mesh * coastline


save = 1
if save == 1:
    dirpath = "./"
    core = "temp_"+f
    outputfile = f"{dirpath}{core}"
    hv.save(p, outputfile, fmt='html', backend='bokeh')

############################make Contour#####################################################################################
#make contor
import numpy as np
## Load an FVCOM model output and plot surface
from PyFVCOM.read import FileReader
from PyFVCOM.read import ncread
from PyFVCOM.plot import Plotter, Time, Depth
from PyFVCOM.tide import make_water_column
from PyFVCOM.grid import unstructured_grid_depths
from cmocean import cm
import matplotlib.pyplot as plt
from pyproj import Proj
import sys
import warnings
# warnings.simplefilter('ignore') # Comment out to show warnings

#fvcom = FileReader(f, dims={'time': slice(0, 1000)}, variables=['zeta', 'temp', 'u'])
#gauge = (140.05, 35.57)  # a sample (lon, lat) position for Estuary example
#index = fvcom.closest_node(gauge).item()  ### ndarray -> scalar
# print(index) # 295

fvcom = FileReader(f, variables=['zeta', 'temp'], dims={'time': slice(0,000), 'node': 150})
time = Time(fvcom, figsize=(40, 9),cmap=cm.thermal, cb_label='{} ({})'.format(fvcom.atts.temp.long_name,
                                                              fvcom.atts.temp.units))
z = unstructured_grid_depths(fvcom.grid.h, fvcom.data.zeta, fvcom.grid.siglay)
## print("z.shape", z.shape)
## fill_seabed makes the part of the plot below the seabed grey.
## We need to squeeze the data array since we've only extracted a single position.
## print(fvcom.grid.h.shape)
## print(fvcom.grid.h)
time.plot_surface(z, fvcom.data.temp, fill_seabed=True, h_offset=0, h_min=fvcom.grid.h)
time.axes.set_ylabel('Depth (m)')
"""






