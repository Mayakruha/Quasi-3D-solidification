import sys
sys.path.append('C:\\Program Files\\ParaView 5.11.1\\bin\\Lib\\site-packages')
#-----------------------------------------------
import numpy as np
from math import exp, tanh, log, pi
import vtk
#------CONSTANTS------------------
AbaqElemTypes={'MASS':1,'CAX4':2,'CPS4':2,'C3D4':3,'C3D6':4,'C3D8':5,'C3D10':6,'C3D15':7,'C3D20R':8,'C3D20':8}
FacesNodes=(None,None,((0,1),(1,2),(2,3),(3,0)),((0,1,2),(0,3,1),(1,3,2),(2,3,0)),((0,1,2),(3,5,4),(0,3,4,1),(1,4,5,2),(2,5,3,0)),\
((0,1,2,3),(4,7,6,5),(0,4,5,1),(1,5,6,2),(2,6,7,3),(3,7,4,0)),\
((0,1,2,4,5,6),(0,3,1,7,8,4),(1,3,2,8,9,5),(2,3,0,9,7,6)),\
((0,1,2,6,7,8),(3,5,4,9,10,11),(0,3,4,1,12,9,13,6),(1,4,5,2,13,10,14,7),(2,5,3,0,14,11,12,8)),\
((0,1,2,3,8,9,10,11),(4,7,6,5,12,13,14,15),(0,4,5,1,16,12,17,8),(1,5,6,2,17,13,18,9),(2,6,7,3,18,14,19,10),(3,7,4,0,19,15,16,11)))
#-----------------------------------------------
#----------------MOULDS-------------------------
#-----------------------------------------------
class mould:
    Height=0.8          #mould height, m    
    Coat_thick0=0.0     #top thickness of coating, m
    Coat_thick1=0.0     #bottom thickness of coating, m
    Tflux_melt=1130     #melting temperature for mould flux, C
    alfa_flux_max=1800  #maximum HTC for flux, W/m2*K
    Rate_flux_drop=0.8  #coefficeint of HTC drop
    Rate_flux_melt=0.05 #coefficient of htan() defining temp range of melting
    Press_coeff=0.022   #coefficient of contact conductence
    Mould_lamda=370     #mould material conductivity, W/mK
    Coat_lamda=80       #coating conductivity, W/mK
    def Prandtl(self,Temp):
        return 12*exp(-0.036*Temp)+1.336343
    # q - heat flux, W
    def HeatUp(self,q):
        self.Twat+=q/self.Qwat*60000/(4230000-3656.2*self.Twat-25.85*self.Twat**2)
    # Twat - Temperature of water, Celsius
    # Qwat - Water flow, l/min
    def Water(self,Twat,Qwat):
        self.Qwat=Qwat
        self.Twat_in=Twat
        self.Twat=self.Twat_in     #current water temeprature
        lamda_wat=0.55748+0.0021525*Twat-0.0000097*Twat**2 #W/m*K
        visc_wat=1.53555258E-06*exp(-0.036*Twat)+2.52805091E-07 #m2/sec
        Pr=self.Prandtl(Twat)
        v_wat=Qwat/60000/self.Wat_sec #m/sec
        Re=v_wat*self.deff/visc_wat
        self.alfa_wat0=0.023*lamda_wat/self.deff*Re**0.8*Pr**0.4*self.geom_coeff(Pr)
        print('\n***Heat transfer parameters:')
        print('Prandtl: '+str(Pr))
        print('Water speed, m/sec: '+str(v_wat))
        print('Reynolds: '+str(Re))
        print('Nominal HTC, kW/(m2K) :'+str(self.alfa_wat0/1000))
        print(' ')
    def geom_coeff(self,Pr):
        return 1
class Circle_tube(mould):
    #---------- Circle ---------------------------
    # D_billet    billet diameter, m
    # Mould_thick distance between water and mould surface, m
    # Wat_thick   thickness of water layer, m    
    def __init__(self,D_billet,Mould_thick,Wat_thick):
        self.R=D_billet/2            # radius of billet cross section, m
        self.Mould_thick=Mould_thick # distance between water and mould surface, m
        self.Wat_thick=Wat_thick     #thickness of water layer, m
        self.deff=2*Wat_thick   # effective diameter, m
        self.Wat_sec=pi*Wat_thick*(2*(D_billet/2+Mould_thick)+Wat_thick)
    def geom_coeff(self, Pr):
        return (1-0.45/(2.4+Pr))*(1+self.Wat_thick/(self.R+self.Mould_thick))**(0.16/Pr**0.15)
    def HeatFlow(self,z,Ts,v,P):#W/m2
        self.alfa_flux=v**0.8*(self.alfa_flux_max*(1-self.Rate_flux_drop/2*(1-tanh(self.Rate_flux_melt*(Ts-self.Tflux_melt))))+self.Press_coeff*P)
        Coat_thick=self.Coat_thick0+(self.Coat_thick1-self.Coat_thick0)*z/self.Height
        self.alfa_wat=self.alfa_wat0*(1+self.Wat_thick/z)
        Q=(Ts-self.Twat)/self.R/(1/self.alfa_wat/(self.R+self.Mould_thick+Coat_thick)+\
            log(1+self.Mould_thick/(self.R+Coat_thick))/self.Mould_lamda+\
            log(1+Coat_thick/self.R)/self.Coat_lamda+\
            1/self.alfa_flux/self.R)
        self.TempW=self.Twat+Q*self.R/(self.R+self.Mould_thick)/self.alfa_wat
        self.alfa_wat*=(self.Prandtl(self.Twat_in)/self.Prandtl(self.TempW))**0.25
        return (Ts-self.Twat)/self.R/(1/self.alfa_wat/(self.R+self.Mould_thick+Coat_thick)+\
                log(1+self.Mould_thick/(self.R+Coat_thick))/self.Mould_lamda+\
                log(1+Coat_thick/self.R)/self.Coat_lamda+\
                1/self.alfa_flux/self.R) #mould 
class Rect_tube(mould):
    #---------- Circle ---------------------------
    # D_billet    billet diameter, m
    # Mould_thick distance between water and mould surface, m
    # Wat_thick   thickness of water layer, m    
    def __init__(self,W_billet,Th_billet,Mould_thick,Wat_thick):
        self.Wdb=W_billet            # width of billet cross section, m
        self.Thb=Th_billet           # thickness of billet cross section, m
        self.Mould_thick=Mould_thick # distance between water and mould surface, m
        self.Wat_thick=Wat_thick     # thickness of water layer, m
        self.deff=2*Wat_thick   # effective diameter, m
        self.Wat_sec=2*Wat_thick*(W_billet+Th_billet+4*Mould_thick+2*Wat_thick)
    def HeatFlow(self,z,Ts,v,P):#W/m2
        self.alfa_flux=v**0.8*(self.alfa_flux_max*(1-self.Rate_flux_drop/2*(1-tanh(self.Rate_flux_melt*(Ts-self.Tflux_melt))))+self.Press_coeff*P)
        Coat_thick=self.Coat_thick0+(self.Coat_thick1-self.Coat_thick0)*z/self.Height
        self.alfa_wat=self.alfa_wat0*(1+self.Wat_thick/z)
        Q=(Ts-self.Twat)/(1/self.alfa_wat+self.Mould_thick/self.Mould_lamda+Coat_thick/self.Coat_lamda+1/self.alfa_flux)
        self.TempW=self.Twat+Q/self.alfa_wat
        self.alfa_wat*=(self.Prandtl(self.Twat_in)/self.Prandtl(self.TempW))**0.25
        return (Ts-self.Twat)/(1/self.alfa_wat+self.Mould_thick/self.Mould_lamda+Coat_thick/self.Coat_lamda+1/self.alfa_flux) #mould 
#-----------------------------------------------
#----------------BILLET-------------------------
#-----------------------------------------------
class Quasi3DTemp:
    Z_end=32.1 #end of calculation, m
    Geom_f=4  #coefficient for section  (1 - full section; 2- half of section etc)
    #---------- casting parameters---------------------------
    v=3.2   #casting speed, m/min
    dTemp=10 #temperature of overheating for liquid steel, K
    Tsr=10   #ambient temperature, Celsius
    MouldLevel=0.1  #mould level,m
    #---------- heat exchange parameters
    Zones=[[0.8,76.74],[1.388,27.54],[3.531,20.82],[5.366,12.54],[8.788,0.3]] #(start of zone, m; water flux, l/m2*min)
    #-----------steel properties----------------------
    lamda=29    #steel conductivity, W/m*K    
#    alfa_liq=2200  #htc for a border between solid and liquid steel, W/m2K
    lamda_liq=190  #conductivity in the area of luquid steel, W/m*K !!!! It doesn't work
    L=272E+3    #heat of solidification, J/kg
    ro_sol=7300 #solid steel density, kg/m3
    ro_liq=7000 #liquid steel density, kg/m3
    Cl=795      #heat capacity for liquid steel, J/kg*K
    Cr=680      #heat capacity for solid steel, J/kg*K
    k=0.2    #power for strain rate in the steel creep equation
    Tc=250      #thermal coefficient in the steel creep equation, K
    Ac=73000    #coefficient for stress in the steel creep equation, MPa
    #-------------------------------------
    wc=(0.125,0.375) #relative coordinates for flux integretaion
    Epsilon=0.001 # tolerance for extraction of temperature from H
    KapaN=20 # Number of columns in the output of billet stiffness
    mould=mould()
    NSets={}
    ESets={}
    Surfs={}  # key-Name of Surface: List (Num of element, Num of face)
    TypeList={}
    MaxNodeNum=0
    MaxElemNum=0
    Coord=np.full(1,None)
    Elems=np.ones(1,dtype=tuple)
    Eltype=np.zeros(1,dtype=np.int8)
    AREA=np.zeros(1)
    T=np.zeros(1)
    NodeLoad={} # key: Name of load; key: Node; Value
    FaceLoad={}
    Faces={} # the first key - min mumber of nodes; the second key - max mumber of nodes
    # list (number of faces, set of numbers of nodes)
    #-----------steel properties----------------------
    def CalcMaterialProp(self,C=0.0,Mn=0.0,Si=0.0,P=0.0,S=0.0,Al=0.0,Cu=0.0,Ni=0.0,Cr=0.0):
        #----------chemical compound of steel, %--------------------------
        self.Tsol=1536-(200*C+16*Si+6*Mn+1.7*Cr+3.9*Ni+93*P+1100*S) #Solidus temperature, Celsius
        self.Tlik=1536-(78*C+7.5*Si+4.9*Mn+1.3*Cr+3.1*Ni+4.7*Cu+3.6*Al+34.4*P+38*S)  #Liquidus temperature, Celsius
        self.Hl=(self.ro_liq*self.Cl+self.ro_sol*self.Cr)*(self.Tlik-self.Tsol)/2+(self.ro_liq+self.ro_sol)*self.L/2 #J/m3
        print('***Thermal Properties of steel:')
        print('Solidus temperature, Celsius: '+str(self.Tsol))
        print('Liquidus temperature, Celsius: '+str(self.Tlik))
        print(' ')
    def Cef(self,Temp):
        if Temp<self.Tsol: return self.Cr
        elif Temp>self.Tlik: return self.Cl
        else: return (self.L+self.Cr*(self.Tlik-Temp)+self.Cl*(Temp-self.Tsol))/(self.Tlik-self.Tsol)
    def FuncTemp(self,Value):
        if Value<self.Tsol: return (Value-self.Tsol)*self.ro_sol*self.Cr
        elif Value>self.Tlik: return (Value-self.Tlik)*self.ro_liq*self.Cl+self.Hl
        else: return (Value-self.Tsol)/(self.Tlik-self.Tsol)*self.Hl
    def Temperature(self,Value):
        if Value<0: return Value/self.ro_sol/self.Cr+self.Tsol
        elif Value>self.Hl: return (Value-self.Hl)/self.ro_liq/self.Cl+self.Tlik
        else: return self.Tsol+(self.Tlik-self.Tsol)*Value/self.Hl
    def Stiffness(self,KapaRate):
        Value=0
        for Node in range(1,self.MaxNodeNum+1):
            if self.T[Node]<self.Tsol:
                Value+=self.Ac*1E+6*(KapaRate*self.Coord[Node][1])**self.k*np.exp(-(self.T[Node]+273)/self.Tc)*self.AREA[Node]*self.Coord[Node][1] # H*m
        return Value/KapaRate # H*m2*sek
    def HeatFlow(self,z,Ts,Name):#W/m2
        if z<=self.Zones[0][0]: return self.mould.HeatFlow(z,Ts,self.v,9.81*self.ro_liq*(z-self.MouldLevel)) #mould 
        else:
            iz=0
            while iz<len(self.Zones)-1 and self.Zones[iz+1][0]<z:iz+=1
            return 142*(Ts-self.Tsr)*self.Zones[iz][1]**0.55*exp(-0.0012*(Ts+273))+5.670367E-8*((Ts+273)**4-(self.Tsr+273)**4) #spray
    def dxydksi(self,nu,xy):
        return -(1-nu)*xy[0]+(1-nu)*xy[1]+nu*xy[2]-nu*xy[3]
    def dxydnu(self,ksi,xy):
        return -(1-ksi)*xy[0]-ksi*xy[1]+ksi*xy[2]+(1-ksi)*xy[3]
#-------------------------------------------------------------------------------
# kj - convergence coefficient (less coefficient, less time step)
# Meshdz - distance between levels for output, [m]
# VTUFileName - name vtu-file for output. There is no output by default
# KapaRate - maximum curvature rate for output of section moment inertia. There is no output by default, [1/m*sec]
# Radius - for visualization, if Radius equals to zero the mesh isn't curved, [m]
    def RunCalc(self, kj=0.5, Meshdz=0.01, VTUFileName='', KapaRate=0, Radius=0):
        face=FacesNodes[self.Eltype[1]][0]
        DistMin=np.linalg.norm(self.Coord[self.Elems[1][face[0]]]-self.Coord[self.Elems[1][face[1]]])
        self.AREA=np.zeros((self.MaxNodeNum+1))
        XEl=np.zeros(4)
        YEl=np.zeros(4)
        for El in range(1,self.MaxElemNum+1):
            for i in range(4):
                XEl[i]=self.Coord[self.Elems[El][i]][0]
                YEl[i]=self.Coord[self.Elems[El][i]][1]
            self.AREA[self.Elems[El][0]]+=(2*XEl[1]*YEl[3]-2*XEl[3]*YEl[1]+(3*XEl[0]-XEl[2])*(YEl[1]-YEl[3])+(3*YEl[0]-YEl[2])*(XEl[3]-XEl[1]))/16
            self.AREA[self.Elems[El][1]]+=(2*XEl[2]*YEl[0]-2*XEl[0]*YEl[2]+(3*XEl[1]-XEl[3])*(YEl[2]-YEl[0])+(3*YEl[1]-YEl[3])*(XEl[0]-XEl[2]))/16
            self.AREA[self.Elems[El][2]]+=(2*XEl[3]*YEl[1]-2*XEl[1]*YEl[3]+(3*XEl[2]-XEl[0])*(YEl[3]-YEl[1])+(3*YEl[2]-YEl[0])*(XEl[1]-XEl[3]))/16
            self.AREA[self.Elems[El][3]]+=(2*XEl[0]*YEl[2]-2*XEl[2]*YEl[0]+(3*XEl[3]-XEl[1])*(YEl[0]-YEl[2])+(3*YEl[3]-YEl[1])*(XEl[2]-XEl[0]))/16
            for face in FacesNodes[self.Eltype[El]]:
                dist=np.linalg.norm(self.Coord[self.Elems[El][face[0]]]-self.Coord[self.Elems[El][face[1]]])
                if DistMin>dist: DistMin=dist
        dtau=kj*DistMin**2*min(self.ro_liq*self.Cl,self.ro_sol*self.Cr)/4/self.lamda  #sek
        dZ=self.v*dtau/60   #m
        Count=int(Meshdz/dZ)
        if Count==0:Count=1
        ElRow=int((self.Z_end-self.MouldLevel)/(dZ*Count))
    #--------------------------------
        if VTUFileName!='':
            Points=vtk.vtkPoints()
            mesh=vtk.vtkUnstructuredGrid()
            mesh.Allocate(ElRow*self.MaxElemNum)
            for i in range(ElRow+1):
                for Node in range(1,self.MaxNodeNum+1):
                    if Radius!=0:
                        if i*dZ*Count<Radius*np.pi/2:
                            Points.InsertNextPoint(self.Coord[Node][0],(self.Coord[Node][1]+Radius)*np.cos(i*dZ*Count/Radius),(self.Coord[Node][1]+Radius)*np.sin(i*dZ*Count/Radius))
                        else:
                            Points.InsertNextPoint(self.Coord[Node][0],Radius*np.pi/2-i*dZ*Count,self.Coord[Node][1]+Radius)
                    else:
                        Points.InsertNextPoint(self.Coord[Node][0],self.Coord[Node][1],i*dZ*Count+self.MouldLevel)
            for i in range(ElRow):
                for El in range(1,self.MaxElemNum+1):
                    mesh.InsertNextCell(vtk.VTK_HEXAHEDRON,8,(self.Elems[El][0]+i*self.MaxNodeNum-1,self.Elems[El][1]+i*self.MaxNodeNum-1,
self.Elems[El][2]+i*self.MaxNodeNum-1,self.Elems[El][3]+i*self.MaxNodeNum-1,
self.Elems[El][0]+(i+1)*self.MaxNodeNum-1,self.Elems[El][1]+(i+1)*self.MaxNodeNum-1,
self.Elems[El][2]+(i+1)*self.MaxNodeNum-1,self.Elems[El][3]+(i+1)*self.MaxNodeNum-1))
            mesh.SetPoints(Points)
            Res=vtk.vtkFloatArray()
            Res.SetName('Temp')
            Res.SetNumberOfValues((ElRow+1)*self.MaxNodeNum)
            Res.Fill(self.Tlik+self.dTemp)
    #--------------------------------
        if KapaRate!=0:
            f=open('BilletStiffness.csv','w')
            for i in range(1,self.KapaN+1):
                f.write(';'+str(i*KapaRate/self.KapaN))
            f.write('\n')
    #--------------------------------
        print('Liquidus Temperature:'+str(self.Tlik)+' C')
        print('Solidus Temperature:'+str(self.Tsol)+' C')
        print('-------------------------------------')
        print('The calculation uses dZ='+str(dZ)+' m')
        print('Number of steps:'+str(ElRow*Count))
        print('-------------------------------------')
        print(str(ElRow+1)+' levels:')
        print('Level: '+str(self.MouldLevel)+' m')
        print('________________________________________________________')
        print(' Axis,m|Twat,C|Twat.wall,C|HTCwat,KW/m2K|Q,kW/m2| Temp,C|')
    #--------------------------------
        ElTemp=np.zeros(4)
        wN=len(self.wc)
        self.T=np.full(self.MaxNodeNum+1,self.Tlik+self.dTemp) #initial temperature, K
        H=np.full(self.MaxNodeNum+1,self.FuncTemp(self.Tlik+self.dTemp))
        Z=self.MouldLevel
        KK=np.zeros((self.MaxElemNum+1,4,4))
        GG=np.zeros((self.MaxElemNum+1,4))
        for El in range(1,self.MaxElemNum+1):
            for i in range(4):
                XEl[i]=self.Coord[self.Elems[El][i]][0]
                YEl[i]=self.Coord[self.Elems[El][i]][1]
            L1=((XEl[0]+XEl[1]-XEl[2]-XEl[3])**2+(YEl[0]+YEl[1]-YEl[2]-YEl[3])**2)**0.5/2
            L2=((XEl[0]+XEl[3]-XEl[1]-XEl[2])**2+(YEl[0]+YEl[3]-YEl[1]-YEl[2])**2)**0.5/2
            GG[El][0]=((XEl[1]-XEl[0])**2+(YEl[1]-YEl[0])**2)**0.5/2 #m
            GG[El][1]=((XEl[2]-XEl[1])**2+(YEl[2]-YEl[1])**2)**0.5/2 #m
            GG[El][2]=((XEl[3]-XEl[2])**2+(YEl[3]-YEl[2])**2)**0.5/2 #m
            GG[El][3]=((XEl[0]-XEl[3])**2+(YEl[0]-YEl[3])**2)**0.5/2 #m
            # n1: ksi=0.5            
            for nu in self.wc:
                KK[El][0]+=dtau*L1/2*np.array((nu-1,1-nu,nu,-nu))/((self.dxydksi(nu,XEl))**2+(self.dxydksi(nu,YEl))**2)**0.5/self.AREA[self.Elems[El][0]]/wN   #sek/m2
                KK[El][1]-=dtau*L1/2*np.array((nu-1,1-nu,nu,-nu))/((self.dxydksi(nu,XEl))**2+(self.dxydksi(nu,YEl))**2)**0.5/self.AREA[self.Elems[El][1]]/wN   #sek/m2
                KK[El][3]+=dtau*L1/2*np.array((nu-0.5,0.5-nu,nu+0.5,-nu-0.5))/((self.dxydksi(nu+0.5,XEl))**2+(self.dxydksi(nu+0.5,YEl))**2)**0.5/self.AREA[self.Elems[El][3]]/wN   #sek/m2
                KK[El][2]-=dtau*L1/2*np.array((nu-0.5,0.5-nu,nu+0.5,-nu-0.5))/((self.dxydksi(nu+0.5,XEl))**2+(self.dxydksi(nu+0.5,YEl))**2)**0.5/self.AREA[self.Elems[El][2]]/wN   #sek/m2
            # n2: nu=0.5            
            for ksi in self.wc:
                KK[El][0]+=dtau*L2/2*np.array((ksi-1,-ksi,ksi,1-ksi))/((self.dxydnu(ksi,XEl))**2+(self.dxydnu(ksi,YEl))**2)**0.5/self.AREA[self.Elems[El][0]]/wN   #sek/m2
                KK[El][3]-=dtau*L2/2*np.array((ksi-1,-ksi,ksi,1-ksi))/((self.dxydnu(ksi,XEl))**2+(self.dxydnu(ksi,YEl))**2)**0.5/self.AREA[self.Elems[El][3]]/wN   #sek/m2
                KK[El][1]+=dtau*L2/2*np.array((ksi-0.5,-ksi-0.5,ksi+0.5,0.5-ksi))/((self.dxydnu(ksi+0.5,XEl))**2+(self.dxydnu(ksi+0.5,YEl))**2)**0.5/self.AREA[self.Elems[El][1]]/wN   #sek/m2
                KK[El][2]-=dtau*L2/2*np.array((ksi-0.5,-ksi-0.5,ksi+0.5,0.5-ksi))/((self.dxydnu(ksi+0.5,XEl))**2+(self.dxydnu(ksi+0.5,YEl))**2)**0.5/self.AREA[self.Elems[El][2]]/wN   #sek/m2
        j=0
        LevelNum=0
        while Z<self.Z_end:
    #----------H calculation---------------
            HeatFlow_level=0
            for El in range(1,self.MaxElemNum+1):
                for i in range(4): ElTemp[i]=self.T[self.Elems[El][i]]
                for i in range(4): H[self.Elems[El][i]]+=np.dot(KK[El][i],ElTemp)*self.lamda  #J/m3
            for BcName in self.Surfs:
                for SurfEl in self.Surfs[BcName]:
                    for El in self.ESets[SurfEl[0]]:                    
                        for i in range(2):
                            Node_indx=FacesNodes[self.Eltype[El]][SurfEl[1]][i]
                            Node=self.Elems[El][Node_indx]
                            Ts=self.T[Node]
                            Q=self.HeatFlow(Z,Ts,BcName)
                            HeatFlow_level+=GG[El][Node_indx]*Q*dZ
                            H[Node]-=dtau*GG[El][Node_indx]*Q/self.AREA[self.Elems[El][Node_indx]] #J/m3
            if Z<=self.mould.Height: self.mould.HeatUp(self.Geom_f*HeatFlow_level)
    #----------Temperature calculation---------------
            for Node in range(1,self.MaxNodeNum+1):self.T[Node]=self.Temperature(H[Node])
            j+=1
            Z+=dZ
            if j==Count:
                print(' {:6.3f}| {:5.2f}|   {:6.2f}  |    {:7.3f}  | {:6.1f}| {:6.1f}|'.format(Z,self.mould.Twat,self.mould.TempW,self.mould.alfa_wat/1000,Q/1000,self.T[1]))
                LevelNum+=1
                if VTUFileName!='':
                    for Node in range(1,self.MaxNodeNum+1):Res.SetValue(self.MaxNodeNum*LevelNum+Node-1,self.T[Node])
                if KapaRate!=0:
                    f.write(str(Z))
                    for i in range(1,self.KapaN+1):
                        f.write(';'+str(self.Stiffness(i*KapaRate/self.KapaN)))
                    f.write('\n')
                j=0          
        if VTUFileName!='':
            mesh.GetPointData().SetScalars(Res)
            output=vtk.vtkXMLUnstructuredGridWriter()
            output.SetInputData(mesh)
            output.SetFileName(VTUFileName)
            output.Write()
        if KapaRate!=0:
            f.close()
#===================================================================
#            READERS:
#===================================================================
#         import Abaqus inp-file
#===================================================================
def import_abq(FileName):
    mesh=Quasi3DTemp()
    mesh.MaxNodeNum=0
    mesh.MaxElemNum=0
    Section=''
    NSet=''
    ESet=''
    ElemNodeNum=0
    ElementType=''
    f=open(FileName,'r')
    txt=f.readline()[:-1]
    while txt:
        if Section=='node':
            while txt and not '*' in txt:
                ValueTxt=txt.split(',')
                NodeNum=int(ValueTxt[0])
                if NodeNum>mesh.MaxNodeNum:mesh.MaxNodeNum=NodeNum
                txt=f.readline()[:-1]
                while '**' in txt: txt=f.readline()[:-1]
            Section=''
        if Section=='element':
            while txt and not '*' in txt:
                ValueTxt=txt.split(',')
                ElemNum=int(ValueTxt[0])
                if ElemNum>mesh.MaxElemNum:mesh.MaxElemNum=ElemNum
                Num=len(ValueTxt)
                if ValueTxt[Num-1]=='':Num-=1                    
                while Num<ElemNodeNum:
                    txt=f.readline()[:-1]
                    ValueTxt=txt.split(',')
                    for Val in ValueTxt:
                        if Val!='':
                            Num+=1                    
                txt=f.readline()[:-1]
                while '**' in txt: txt=f.readline()[:-1]
            Section=''
        if '*node' in txt.lower() and not '*node output' in txt.lower(): Section='node'
        if '*element' in txt.lower() and not '*element output' in txt.lower():
            Section='element'
            SetNamePos=txt.lower().find('type')+5
            if ',' in txt[SetNamePos:]: ElementType=txt[SetNamePos:txt.find(',',SetNamePos)]
            else: ElementType=txt[SetNamePos:]
            mesh.TypeList[AbaqElemTypes[ElementType]]=ElementType
            if AbaqElemTypes[ElementType]!=None:
                ElemNodeNum=max(max(FacesNodes[AbaqElemTypes[ElementType]]))+1
            else:ElemNodeNum=1
        txt=f.readline()[:-1]
    f.close()
    mesh.Coord=np.full(mesh.MaxNodeNum+1,None)
    mesh.Elems=np.ones(mesh.MaxElemNum+1,dtype=tuple)
    mesh.Eltype=np.zeros(mesh.MaxElemNum+1,dtype=np.int8)
    f=open(FileName,'r')
    Section=''
    NSet=''
    ESet=''
    Surf=''        
    txt=f.readline()[:-1]
    while txt:
        if Section=='node':
            while txt and not '*' in txt:
                ValueTxt=txt.split(',')
                NodeNum=int(ValueTxt[0])
                if NSet!='': mesh.NSets[NSet].append(NodeNum)
                mesh.Coord[NodeNum]=np.array(list(map(float,ValueTxt[1:])))
                txt=f.readline()[:-1]
                while '**' in txt: txt=f.readline()[:-1]
            Section=''
            NSet=''
        elif NSet!='':
            while txt and not '*' in txt:
                for Val in txt.replace(' ','').split(','):
                    if Val in mesh.NSets:
                        for NodeNum in mesh.NSets[Val]: mesh.NSets[NSet].append(NodeNum)
                    elif Val!='': mesh.NSets[NSet].append(int(Val))
                txt=f.readline()[:-1]
                while '**' in txt: txt=f.readline()[:-1]
            NSet=''           
        if Section=='element':
            while txt and not '*' in txt:
                ValueTxt=txt.split(',')
                ElemNum=int(ValueTxt[0])
                mesh.Eltype[ElemNum]=AbaqElemTypes[ElementType]
                if ESet!='': mesh.ESets[ESet].append(ElemNum)
                Num=len(ValueTxt)
                if ValueTxt[Num-1]=='':Num-=1
                mesh.Elems[ElemNum]=list(map(int,ValueTxt[1:Num]))
                while Num<ElemNodeNum:
                    txt=f.readline()[:-1]
                    ValueTxt=txt.split(',')
                    for Val in ValueTxt:
                        if Val!='':
                            mesh.Elems[ElemNum].append(int(Val))
                            Num+=1
                txt=f.readline()[:-1]
                while '**' in txt: txt=f.readline()[:-1]
            Section=''
            ESet=''
            ElementType=''
        elif ESet!='':
            while txt and not '*' in txt:
                for Val in txt.replace(' ','').split(','):
                    if Val in mesh.ESets:
                        for ElemNum in mesh.ESets[Val]: mesh.ESets[ESet].append(ElemNum)
                    elif Val!='': mesh.ESets[ESet].append(int(Val))                   
                txt=f.readline()[:-1]
                while '**' in txt: txt=f.readline()[:-1]
            ESet=''
        if Surf!='':
            while txt and not '*' in txt:
                ValueTxt=txt.replace(' ','').split(',')
                mesh.Surfs[Surf].append((ValueTxt[0],int(ValueTxt[1][1:])-1))               
                txt=f.readline()[:-1]
                while '**' in txt: txt=f.readline()[:-1]
            Surf=''            
        if '*node' in txt.lower() and not '*node output' in txt.lower():
            Section='node'
            txt.replace(' ','')
            SetNamePos=txt.lower().find('nset')
            if SetNamePos>4:
                SetNamePos=txt.find('=',SetNamePos)
                if ',' in txt[SetNamePos:]: NSet=txt[SetNamePos+1:txt.find(',',SetNamePos)]
                else: NSet=txt[SetNamePos+1:]                    
                if not NSet in mesh.NSets: mesh.NSets[NSet]=[]
        if '*element' in txt.lower() and not '*element output' in txt.lower():
            Section='element'
            txt.replace(' ','')
            SetNamePos=txt.lower().find('elset')
            if SetNamePos>4:
                SetNamePos=txt.find('=',SetNamePos)
                if ',' in txt[SetNamePos:]: ESet=txt[SetNamePos+1:txt.find(',',SetNamePos)]
                else: ESet=txt[SetNamePos+1:]
                if not ESet in mesh.ESets: mesh.ESets[ESet]=[]
            SetNamePos=txt.lower().find('type')+5
            if ',' in txt[SetNamePos:]: ElementType=txt[SetNamePos:txt.find(',',SetNamePos)]
            else: ElementType=txt[SetNamePos:]
            if AbaqElemTypes[ElementType]!=None:
                ElemNodeNum=max(max(FacesNodes[AbaqElemTypes[ElementType]]))+1
            else:ElemNodeNum=1
        if '*nset' in txt.lower():
            txt.replace(' ','')
            SetNamePos=txt.lower().find('nset',3)
            SetNamePos=txt.find('=',SetNamePos)
            if ',' in txt[SetNamePos:]: NSet=txt[SetNamePos+1:txt.find(',',SetNamePos)]
            else: NSet=txt[SetNamePos+1:]
            if not NSet in mesh.NSets: mesh.NSets[NSet]=[]
            if 'generate' in txt.lower():
                txt=f.readline()[:-1]
                ValueTxt=txt.split(',')
                for NodeNum in range(int(ValueTxt[0]),int(ValueTxt[1])+int(ValueTxt[2]),int(ValueTxt[2])):
                    mesh.NSets[NSet].append(NodeNum)
                NSet=''
        if '*elset' in txt.lower():
            txt.replace(' ','')
            SetNamePos=txt.lower().find('elset',4)
            SetNamePos=txt.find('=',SetNamePos)
            if ',' in txt[SetNamePos:]: ESet=txt[SetNamePos+1:txt.find(',',SetNamePos)]
            else: ESet=txt[SetNamePos+1:]
            if not ESet in mesh.ESets: mesh.ESets[ESet]=[]
            if 'generate' in txt.lower():
                txt=f.readline()[:-1]
                ValueTxt=txt.split(',')
                for ElemNum in range(int(ValueTxt[0]),int(ValueTxt[1])+int(ValueTxt[2]),int(ValueTxt[2])):
                    mesh.ESets[ESet].append(ElemNum)
                ESet=''
        if '*surface' in txt.lower():
            txt.replace(' ','')
            SetNamePos=txt.lower().find('name',7)
            if ',' in txt[SetNamePos:]: Surf=txt[SetNamePos+5:txt.find(',',SetNamePos)]
            else: Surf=txt[SetNamePos+5:]
            if not Surf in mesh.Surfs: mesh.Surfs[Surf]=[]
        txt=f.readline()[:-1]
    f.close()
    print('Model has been imported')
    return mesh
