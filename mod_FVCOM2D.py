import netCDF4
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import holoviews as hv
import datashader.utils as du

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
        #self.triang = tri.Triangulation(self.x, self.y, triangles=nv)
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
        #self.nbe = np.array([[self.nv[n, j], self.nv[n, (j+2)%3]] \
        #    for n in range(len(self.triang.neighbors)) \
        #    for j in range(3) if self.triang.neighbors[n,j] == -1])

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
