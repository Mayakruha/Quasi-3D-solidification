import sys
sys.path.append('D:\\Scripts\\FEMtoolkit\\')
sys.path.append('C:\\Program Files\\ParaView 5.11.1\\bin\\Lib\\site-packages\\')
from FEMtoolkit import import_abq
#Sizes={156:160,210:217,228:235,260:268,310:320,330:340,340:367,350:361,360:371,410:423}
#for Size in Sizes:
#    model=import_abq('D://VTP//CCM//3DCCM//Sec_196.inp')
#    model.Scale(Sizes[Size]/206)
#    model.export_abq('D://VTP//CCM//3DCCM//Sec_'+str(Size)+'.inp')
#---------------------------------------------------------------
model=import_abq('D://VTP//CCM//3DCCM//Sec_360x360.inp')
model.NodeIntoSurf('Outer')
mesh=model.ReNumb()
mesh.NSets={}
mesh.export_abq('D://VTP//CCM//3DCCM//Sec_360x360-1.inp')