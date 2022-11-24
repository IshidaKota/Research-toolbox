#for ITO
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
###

import netCDF4
import matplotlib.pyplot as plt
import sys
import numpy as np


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
stn = ['chiba1buoy','chibaharo','kawasaki','urayasu']
stn_node = [136,46,222,51]
ncfile = sys.argv[1]
path = '../run/'+ ncfile+'_0001.nc'
nc = netCDF4.Dataset(path,'r')
temp = nc.variables['temp']
"""
j = 0
for index in stn_node:
    fig = plt.figure(figsize=(10,4))
    print('start plotting')
    ax = fig.add_subplot(1,1,1,xlabel = 'days since 2020/1/1 0:00:00(UTC)',ylabel='temprature degree')
    for i in range(0,19,5):
        ax.plot(temp[:,i, index],linewidth=0.5,label='sigma{}'.format(i)) #float64 temp(time, siglay, node)
    fig.suptitle('temprature@{}'.format(stn[j]))
    ax.legend()
    plt.grid()
    plt.ylim(0,40)
    plt.xticks([24*30*i for i in range(12)], [30*i for i in range(12)])
    plt.plot()
    fig.savefig('./png/temp_'+sys.argv[1]+'_'+stn[j]+'.png')
    print('done,{}'.format(stn[j]))
    j+=1
    
"""



#make contor

t =nc.variables['time'][:]


for index in stn_node:
    sig = [1,2]
    #sig = [nc.variables['zeta'][index]*i/20 for i in range(19)]
    T,S = np.meshgrid(t,sig)
    Temp = np.zeros((len(t),len(sig)))
    for i in range(len(t)):
        for j in range(len(sig)):
            Temp[i,j] = temp[i,j,index]
    #cont=plt.contour(T,S,Temp,  5, vmin=5,vmax=20, colors=['black'])
    #cont.clabel(fmt='%1.1f', fontsize=14)
    def tempr(T,S,Temp):
        return Temp
    

    #plt.xlabel('time', fontsize=24)
    #plt.ylabel('siglay', fontsize=24)


    plt.pcolormesh(T,S,tempr(T,S,Temp), shading='flat',cmap='rainbow') #カラー等高線図
    pp=plt.colorbar (orientation="vertical") # カラーバーの表示 
    pp.set_label("Label",  fontsize=24)

    plt.plot()
    fig.savefig('./png/temp_contour'+str(index)+'.png')

nc.close()
