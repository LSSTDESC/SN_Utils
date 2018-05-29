import numpy as np

class Observations:
    def __init__(self,ra=-1,dec=-1,filename='',nseasons=10,season_length=95,names=['mjd','Ra','Dec']):
        
        self.mjd=names[0]
        self.ra=names[1]
        self.dec=names[2]

        if filename != '':
            data=self.Load(filename)
            #print(data.dtype)
            #data=np.sort(data,order='mjd')
            data.sort(order=self.mjd)
            self.all_seasons=data
            self.seasons=self.Get_Seasons(data,season_length)
            self.Remove_Poor_Seasons()
            self.Ra=np.unique(data[self.ra])[0]
            self.Dec=np.unique(data[self.dec])[0]
        else:
            self.Ra=ra
            self.Dec=dec
            self.seasons={}
            for i in range(nseasons):
                self.seasons[i]=None

        #print 'Nseasons',len(self.seasons)
        #print data

    def Load(self,filename):

        sfile=open(filename,'r')
        varname=[]
        r=[]
        for line in sfile.readlines():
            if line[0] == '#':
                varname.append(line.split(' ')[1])
            else:
                tofill=[]
                thesplit=line.strip().split(' ')
                tofill.append(thesplit[0])
                for i in range(1,len(thesplit)):
                    tofill.append(float(thesplit[i]))
                r.append(tuple(tofill))
        return np.rec.fromrecords(r,names=varname)
            
    def Get_Seasons(self,data,season_length):
        
        thediff=data[self.mjd][1:]-data[self.mjd][:-1]
        idx,=np.where(thediff > season_length)
        lidx=[val+1 for val in list(idx)]
        
        lidx.insert(0,0)
        lidx.append(len(data[self.mjd]))
        #print 'ay',lidx
        
        seasons={}
        
        for i in range(len(lidx)-1):
            #print lidx[i],lidx[i+1]
            seasons[i]=data[lidx[i]:lidx[i+1]]

        return seasons
            

    def Remove_Poor_Seasons(self):

        iseason=-1
        res={}

        for i in range(len(self.seasons)):
            #print 'season',i,len(self.seasons[i])
            if len(self.seasons[i])>10:
                iseason+=1
                res[iseason]=self.seasons[i]

        self.seasons=res



    
