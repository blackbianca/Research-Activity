import matplotlib.pyplot as p
import numpy as n
import matplotlib as m
from scipy import optimize
m.rcParams.update({'font.size':7})


#polinomyal function to fit the data
def pol(x,a0,a1,a2,a3):
    #y = a0 + a1*x + a2*x**2. + a3*x**3. + a4*x**0.5 + a5*x**1.5
    y = a0 + a1*x + a2*x**2. + a3*x**3. 
    return y




def read():
        f606, f814, flag,rms1,rms2=n.genfromtxt('horologium i.txt', dtype=float, comments='#', usecols=(2,5,10,3,6), unpack=True)
        f6=[]
        f8=[]
        i=0
        for i in range(len(f606)):
                if (rms1[i]<=1):
                        if (rms2[i]<=1):
                                if (flag[i]==1):
                                        f6.append(f606[i])
                                        f8.append(f814[i])
                i+=1
        return f6,f8

f6,f8=read()
f606=n.array(f6)
f814=n.array(f8)
N=len(f606)
color=n.zeros(N,float)

i=0
for i in range(len(color)):
        color[i]=f606[i]-f814[i]
        i+=1
        

########################## iso 13.5 alpha 0.8#################

f6_13,f8_13=n.genfromtxt('13,5-0.8.iso', dtype=float, comments='#', usecols=(10,15), unpack=True)
N5=int(len(f6_13))
color_iso13=n.zeros(N5,float)
w=0
for w in range(N5):
        color_iso13[w]=f6_13[w]-f8_13[w]
        w+=1
        
###f8 filter        
c_excess8_13=0.0146    
d8_13=19.7
A8_13=1.842*c_excess8_13
ident8_13=n.ones(N5,float)
M813=n.zeros(N5,float)
M813[:]=f8_13[:]+d8_13*ident8_13[:]+A8_13*ident8_13[:]

### f6 filter
c_eccess6_13=0.0146
d6_13=19.7
A6_13=2.8782*c_eccess6_13
ident6_13=n.ones(N5,float)
M613=n.zeros(N5,float)
M613[:]=f6_13[:]+d6_13*ident6_13[:]+A6_13*ident6_13[:]





#################iso 14   alpha 0.4         fe/H=-2.49 ###########################
f6_14,f8_14=n.genfromtxt('14Gyrs.iso', dtype=float, comments='#', usecols=(10,15), unpack=True)
N3=int(len(f6_14))
color_iso14=n.zeros(N3,float)

for j in range(N3):
        color_iso14[j]=f6_14[j]-f8_14[j]
        j+=1
        
###f8 filter
c_excess8_14=0.0146    
d8_14=19.7
A8_14=1.842*c_excess8_14
ident8_14=n.ones(N3,float)
M814=n.zeros(N3,float)
M814[:]=f8_14[:]+d8_14*ident8_14[:]+A8_14*ident8_14[:]

###f6 filter
d6_14=19.7
c_excess6_14=0.0146
A6_14=2.8782*c_excess6_14
ident6_14=n.ones(N3,float)
M614=n.zeros(N3,float)
M614[:]=f6_14[:]+d6_14*ident6_14[:]+A6_14*ident6_14[:]



########################## iso 13.7 alpha 0.2   fe=-2.4 #################

f6_137,f8_137=n.genfromtxt('13,7-0,2.iso', dtype=float, comments='#', usecols=(10,15), unpack=True)
N6=int(len(f6_137))
color_iso137=n.zeros(N6,float)
l=0
for l in range(N6):
        color_iso137[l]=f6_137[l]-f8_137[l]
        l+=1

###f8 filter        
c_excess8_137=0.0146    
d8_137=19.7
A8_137=1.842*c_excess8_137
ident8_137=n.ones(N6,float)
M8137=n.zeros(N6,float)
M8137[:]=f8_137[:]+d8_137*ident8_137[:]+A8_137*ident8_137[:]

### f6 filter
c_eccess6_137=0.0146
d6_137=19.7
A6_137=2.8782*c_eccess6_137
ident6_137=n.ones(N6,float)
M6137=n.zeros(N6,float)
M6137[:]=f6_137[:]+d6_137*ident6_137[:]+A6_137*ident6_137[:]


############################### PLOTS ###################
#plot f8 filter
p.plot(color,f814, 'b.', label='data')
p.plot(color_iso13,M813,color='red', label='iso 13.5 Gyrs')
p.plot(color_iso14,M814, color='yellow', label='iso 14 Gyrs')
p.plot(color_iso137,M8137,color='aqua', label='iso 13.7 Gyrs')
p.legend()
p.xlabel('m606-m814')
p.ylabel('m814')
p.axis([0,1.4,26,18])
p.show()

##plot f6 filter
p.plot(color,f606, 'c.', label='data')
p.plot(color_iso137,M6137,color='magenta',label='iso 13.7 Gyrs')
p.plot(color_iso14,M614,'b-', label='iso 14 Gyrs')
p.plot(color_iso13,M613,color='red',label='iso 13.5 Gyrs')
p.legend()
p.axis([0,1.4,26,18])
p.xlabel('m606-m814')
p.ylabel('m606')
p.show()

############################### BINARIES ########################

mass1,m16,m18=n.genfromtxt('13,7-0,2.iso', dtype=float, comments='#', usecols=(1,10,15), unpack=True)

M16=[]
M26=[]
M18=[]
M28=[]
#### l j w i 
e=1e-3
for k in range(0,245):
        mass16=mass1[k]
        mass26=mass16/2
        for t in range(len(mass1)):
                if (abs(mass1[t]-mass26)<e):
                        M16.append(m16[k])
                        M26.append(m16[t])
                        M18.append(m18[k])
                        M28.append(m18[t])
        
#for q in range(len(mass1)):
 #       mass18=mass1[q]
  #      mass28=mass18/2
   #     for o in range(len(mass1)):
    #            if (abs(mass1[o]-mass28)<e):
                        

m16=n.array(M16)
m18=n.array(M18)
m26=n.array(M26)
m28=n.array(M28)                      
N8=int(len(m16))
N9=int(len(m18))

mag6=n.zeros(N8,float)
mag6[:]=-2.5*n.log10((10**(-m16[:]/2.5))+10**(-m26[:]/2.5))
mag8=n.zeros(N9,float)
mag8[:]=-2.5*n.log10((10**(-m18[:]/2.5))+10**(-m28[:]/2.5))       

M6b=n.zeros(N8,float)
M8b=n.zeros(N9, float)

identt6=n.ones(N8,float)
identt8=n.ones(N9,float)

M6b[:]= mag6[:]+d6_137*identt6[:]+A6_137*identt6[:]
M8b[:]=mag8[:]+d8_137*identt8[:]+A8_137*identt8[:]

       
col=n.zeros(N9, float)
col[:]=mag6[:]-mag8[:]

################## interpolation main sequence

ms6 = []
ms8 = []
colms6 = []
colms8 = []

for i in range(int(len(M6137))):
    if(M6137[i]>23.5 and M6137[i]<30):
        ms6.append(M6137[i])
        colms6.append(color_iso137[i])
for i in range(int(len(M8137))):
    if(M8137[i]>23.5 and M8137[i]<30):
        ms8.append(M8137[i])
        colms8.append(color_iso137[i])


par = [1.,1.,1.,1.,1.,1.,1.]
poptms, pcovms = optimize.curve_fit(pol, colms6, ms6,p0=(par[0], par[1], par[2], par[3]))
popt2ms, pcov2ms = optimize.curve_fit(pol, colms8, ms8,p0=(par[0], par[1], par[2], par[3]))

ctest = n.linspace(0.4, 1., int(1e4))
ms6test = pol(ctest, poptms[0], poptms[1], poptms[2], poptms[3])
ms8test = pol(ctest, popt2ms[0], popt2ms[1], popt2ms[2], popt2ms[3])

################################# line q=1

M6bint1 = ms6test - 0.75
M8bint1 = ms8test - 0.75

################ linear interpolation of q=0.5
M6btest=[]
M8btest=[]
ctest6=[]
ctest8=[]
for i in range(int(len(M6b))):
        if (M6b[i]>23.5):
                M6btest.append(M6b[i])
                ctest6.append(col[i])
for i in range(int(len(M8b))):
        if (M8b[i]>23.5):
                M8btest.append(M8b[i])
                ctest8.append(col[i])


par = [1.,1.,1.,1.,1.]
popt, pcov = optimize.curve_fit(pol, ctest6, M6btest,p0=(par[0], par[1], par[2], par[3]))
popt2, pcov2 = optimize.curve_fit(pol, ctest8, M8btest,p0=(par[0], par[1], par[2], par[3]))

M6bint=pol(ctest, popt[0], popt[1], popt[2], popt[3])
M8bint=pol(ctest, popt2[0], popt2[1], popt2[2], popt2[3])


######################## filling regions in filter 606
regAx = []
regAy = []
regBx = []
regBy = []

for i in range(int(len(f6))):
        if (f6[i]>23.5 and f6[i]<30.):
                c = color[i]
                ms = pol(c, poptms[0], poptms[1], poptms[2], poptms[3])
                q1 = ms-0.75
                if(f6[i]<ms and f6[i]>q1):
                #if(f6[i]>q1):
                        regAx.append(c)
                        regAy.append(f6[i])
                        q05 = pol(c, popt[0], popt[1], popt[2], popt[3])
                        if(f6[i]<q05):
                                regBx.append(c)
                                regBy.append(f6[i])

freq6 = float(int(len(regBx))/int(len(regAx)))
print(freq6)

                        
################# PLOTS ####################
p.plot(color,f606, 'c.', label='data')
p.plot(color_iso137,M6137,color='magenta',label='iso 13.7 Gyrs')
#p.plot(color_iso14,M614,'b-', label='iso 14 Gyrs')
#p.plot(color_iso13,M613,color='red',label='iso 13.5 Gyrs')
#p.plot(ctest,m6bintest,color='yellow',label='bin q=0.5')
#p.plot(color_iso137,Mq16,color='red',label='bin q=1')
p.scatter(regAx, regAy, color="orange", label="region A")
p.scatter(regBx, regBy, color="green", label="region B")
p.plot(ctest,ms6test,color='blue',label='Main Sequence fit')
p.plot(ctest,M6bint1, 'black', label="q=1 line")
p.plot(col,M6b, 'y-', label="bin q=0.5")
p.plot(ctest,M6bint,color='red',label='q=0.5 line')
p.legend()
p.axis([0,1.4,26,18])
p.xlabel('m606-m814')
p.ylabel('m606')

######################## filling regions in filter 814
regAx8 = []
regAy8 = []
regBx8 = []
regBy8 = []

for i in range(int(len(f8))):
        if (f8[i]>23.5 and f8[i]<30.):
                c = color[i]
                ms = pol(c, popt2ms[0], popt2ms[1], popt2ms[2], popt2ms[3])
                q1 = ms-0.75
                if(f8[i]<ms and f8[i]>q1):
                #if(f8[i]>q1):
                        regAx8.append(c)
                        regAy8.append(f8[i])
                        q05 = pol(c, popt2[0], popt2[1], popt2[2], popt2[3])
                        if(f8[i]<q05):
                                regBx8.append(c)
                                regBy8.append(f8[i])

freq8 = float(int(len(regBx8))/int(len(regAx8)))
print(freq8)


p.show()

p.plot(color,f814, 'c.', label='data')
#p.plot(color_iso13,M813,color='red', label='iso 13.5 Gyrs')
#p.plot(color_iso14,M814, color='yellow', label='iso 14 Gyrs')
p.plot(color_iso137,M8137,color='aqua', label='iso 13.7 Gyrs')
#p.plot(ctest,m8bintest,color='yellow',label='bin q=0.5')
#p.plot(color_iso137,Mq18,color='red', label='bin q=1')
p.scatter(regAx8, regAy8, color="orange", label="region A")
p.scatter(regBx8, regBy8, color="green", label="region B")
p.plot(ctest,ms8test,color='blue',label='Main Sequence fit')
p.plot(ctest,M8bint1, 'black', label="q=1 line")
p.plot(col,M8b, 'y-', label="bin q=0.5")
p.plot(ctest,M8bint,color='red',label='q=0.5 line')
p.legend()
p.xlabel('m606-m814')
p.ylabel('m814')
p.axis([0,1.4,26,18])
p.show()














