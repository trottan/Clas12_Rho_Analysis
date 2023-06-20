
h_p = [] 

tCount = 1
QCount = 1
phiCount = 0
for i in phiBins:
    if not (QCount == 1 and (tCount == 2 or tCount == 3) and (phiCount == 0 or phiCount == 8)):
        h_p.append(i.Filter("ihel > 0").Histo1D(("hmrho","IM:#pi^{+}+#pi^{-},Heli+, Qbin:" + str(QCount) + "tBin" + str(tCount) + "phibin:" + str(phiCount),200,0,2.5), "mrho"))
    #print(str(tCount) + "phibin:" + str(phiCount))
    phiCount +=1 
    if phiCount % 9 ==0:
        phiCount = 0
        tCount += 1
    if tCount == 4:
        tCount = 1
        QCount +=1

#hp,rhoAmp_p, rhoAmpError_p = fit1d(h_p)



h_m = [] 

QCount = 1
tCount = 1
phiCount = 0
for i in phiBins:
    if not (QCount == 1 and (tCount == 2 or tCount == 3) and (phiCount == 0 or phiCount == 8)):
        h_m.append(i.Filter("ihel < 0").Histo1D(("hmrho","IM:#pi^{+}+#pi^{-},Heli- Qbin:" + str(QCount) + "tBin" + str(tCount) + "phibin:" + str(phiCount),200,0,2.5), "mrho"))
    phiCount +=1 
    if phiCount % 9 ==0:
        phiCount = 0
        tCount += 1
    if tCount == 4:
        tCount = 1
        QCount +=1

#hm,rhoAmp_m, rhoAmpError_m = fit1d(h_m)

print("len is",len(h_p),len(h_m))

h = h_p + h_m


hc,rhoAmp, rhoAmpError = fit1d(h)

#print(rhoAmp)
#print()
#print(rhoAmpError)


rhoAmp_p,rhoAmp_m = np.split(np.array(rhoAmp),2)
rhoAmpError_p,rhoAmpError_m = np.split(np.array(rhoAmpError),2)

print(repr(rhoAmp_p))
print(repr(rhoAmp_m))
print()
print(repr(rhoAmpError_p))
print(repr(rhoAmpError_m))
