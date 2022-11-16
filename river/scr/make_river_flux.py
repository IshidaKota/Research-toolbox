#author ishid 2022/11/15


from glob import glob
import pandas as pd
from datetime import datetime as dt
from datetime import timedelta
import codecs

def make_river_flux(place_list):
    """
水質水文データベースから各月ごとにdat形式でダウウンロードしたデータをfvcom-toolbox
で処理できる形のcsvに出力する。datファイル名に統一した名前を入れ、その名前をlistで
引数に入れる。
    """
    for place in place_list:
        paths =glob('../'+place+'*') #datがあるpathを指定してください！
        
        #一日のデータをresの各行で格納
        res = []
        for path in paths:
            with codecs.open(path,'r','utf-8','ignore') as f:
                l = f.readlines(encoding='cp932')
                for i in range(10,len(l)): #11行目以降が時系列データになっています。
                    res.append(l[i].split(','))
        


        #一日のデータのうち、空白や記号を取り除きます。奇数番目にデータが入っている
        data = []
        for i in range(len(res)):
            data_cld = [res[i][j].replace('    ','') for j in range(1,len(res[i]),2)]
            data.append(data_cld)
        
        #一行のlistにまとめる
        flux = [float(data[i][j]) for i in range(len(res))for j in range(24)]

        datelist = []
        for x in range((len(res))*24):
            t = dt.strptime(res[0][0], '%Y/%m/%d')+ timedelta(hours=1) \
                    +timedelta(hours=x) #初期時間からhours=xだけずらした時間
            t = t.strftime('%d-%m-%Y %H:%M') #書式の変換
            datelist.append(t)
            
        #create df
        df = pd.DataFrame()
        df['datetime'] = datelist
        df[place] = flux
        df.to_csv('../river_flux2020.csv')
        return df

flux = make_river_flux(['tamagawa'])