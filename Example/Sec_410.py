from Quasi3DTemp import import_abq
model=import_abq('Sec_410.inp')
model.v=0.45 #m/min
model.mould=Circle_tube(0.423,0.026,0.0035) # billet diameter [m], distance between water and mould surface [m], thickness of water layer [m]
model.mould.Water(10,1900) #Temperature of water [Celsius], Water flow [l/min]
model.Zones=[[0.8,85.8],[1.388,30.78],[3.531,23.28],[5.366,13.98],[8.788,0.3]] #(start of zone [m]; water flux [l/m2*min])
model.CalcMaterialProp(C=1.05,Mn=1.2,Si=0.65,P=0.027,S=0.02,Cu=0.25,Ni=0.3,Cr=1.65)
model.RunCalc(VTUFileName='Temp_410.vtu', Radius=14, Meshdz=0.1, KapaRate=0.0005)
