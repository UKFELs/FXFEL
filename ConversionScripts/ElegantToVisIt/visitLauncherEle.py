def SetupSerialLauncher():
  lp2=visit.LaunchProfile()
  lp2.SetTimeout(60)
  lp2.SetTimeLimitSet(1)
  lp2.SetTimeLimit("60")
  lp2.SetProfileName("serial")
  lp2.SetActive(1)
  lp2.SetArguments("-debug 5")
  lp2.SetLaunchMethod('serial')
  return lp2

def SetupPhase2():
  phase2=visit.MachineProfile()
  phase2.SetTunnelSSH(1)
  phase2.SetUserName('jds89-bwm06')
  phase2.SetDirectory('/gpfs/stfc/local/HCP084/bwm06/shared/visit2_10_0.linux-x86_64')
  phase2.SetActiveProfile(1)
  phase2.SetHost('phase2.wonder.hartree.stfc.ac.uk')
  phase2.SetClientHostDetermination(2)
  lp2=SetupSerialLauncher()
  phase2.AddLaunchProfiles(lp2)
  phase2.SetActiveProfile(0)
  ce=visit.OpenComputeEngine(phase2)
  return phase2

# either load parameters from file, or enter them afresh using serialize/pickle
# LocalSettings
localVisItDir= "/home/jonny/visit2_10_0.linux-x86_64/2.10.0"
localPythonPackageDir = "/home/jonny/visit2_10_0.linux-x86_64/2.10.0/linux-x86_64/lib/site-packages" 
runRemotely=0
remoteTimeSeriesAstraDB='phase2.wonder.hartree.stfc.ac.uk:/gpfs/stfc/local/HCP084/bwm06/shared/testastra/working/test_analyzeBeam.h5'
remoteTimeSeriesEleSigmaDB='phase2.wonder.hartree.stfc.ac.uk:/gpfs/stfc/local/HCP084/bwm06/shared/testele/phase100/clara.sig-page1.h5'
remoteTimeSeriesElePhaseSpaceDB='phase2.wonder.hartree.stfc.ac.uk:/gpfs/stfc/local/HCP084/bwm06/shared/testele/clara_track_electrons_*.vsh5'
#Remote settings
def isOpenDatabase(filename):
  pass
  return 0

def PhaseSpace(x1,x2,dumpNo):
   ScatAttrs = visit.ScatterAttributes()
   ScatAtts.var1 = x1
   ScatAtts.var2 = x2
   ScatAtts.var1Role = ScatterAtts.Coordinate0
   ScatAtts.var2Role = ScatterAtts.Coordinate1
   ScatAtts.var3Role = ScatterAtts.None
   ScatAtts.var4Role = ScatterAtts.None
   visit.AddPlot('Scatter',x1,1,1) 

def PhaseSpaceXY(elementNo=1):
  if isOpenDataBase(''):
    PhaseSpace('electrons_x','electrons_y')
  else:
    psData=visit.OpenDatabase(,elementNo,'Vs')
    PhaseSpace('electrons_x','electrons_y')
    
import sys
sys.path.insert(0,localPythonPackageDir)
import visit
visit.Launch(vdir=localVisItDir)
if runRemotely:
  p2=SetupPhase2()
  data2=visit.OpenDatabase(remoteTimeSeriesAstraDB,0,'Vs')
else:
  data=visit.OpenDatabase
visit.AddPlot('Curve','s1')
visit.AddPlot('Curve','s2')
visit.DrawPlots()
visit.AddWindow()
visit.AddPlot('Curve','pAverage')
visit.DrawPlots()
visit.SetActiveWindow(1)
AnnotationAtts = visit.AnnotationAttributes()
AnnotationAtts.axes2D.yAxis.grid = 1
AnnotationAtts.axes2D.xAxis.grid = 1
AnnotationAtts.axes2D.xAxis.title.userTitle = 1
AnnotationAtts.axes2D.xAxis.title.userUnits = 1
AnnotationAtts.axes2D.xAxis.title.title = "Position along the machine s"
AnnotationAtts.axes2D.xAxis.title.units = "m"
visit.SetAnnotationAttributes(AnnotationAtts)
visit.OpenGUI()
