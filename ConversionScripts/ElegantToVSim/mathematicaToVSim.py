import csv
fh=open("../Manchester/ElectronPlasmaWake/PARS_Beams/VorpalTestData_20pC_1000_Elegant.dat","rb")
myreader=csv.reader(fh,delimiter="\t")
fhout=open("parsUltraShort20pC10002d.dat","wb")
mywriter=csv.writer(fhout,delimiter=" ")
c=299792458 #m/s
for thisrecord in myreader:
  VorpalX=-float(thisrecord[2])*c
  VorpalPX=float(thisrecord[5])*c
  VorpalY=float(thisrecord[1])
  VorpalZ=float(thisrecord[0]) # Throw this away for a 2D sim.
  VorpalPY=float(thisrecord[4])*VorpalPX
  VorpalPZ=float(thisrecord[3])*VorpalPX
  mywriter.writerow([VorpalX,VorpalY,VorpalPX,VorpalPY,VorpalPZ])
fh.close()
fhout.close()
