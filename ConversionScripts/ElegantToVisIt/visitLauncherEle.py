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

def SaveWindow(filename):
  SaveWindowAtts = visit.SaveWindowAttributes()
  SaveWindowAtts.outputToCurrentDirectory = 1
  SaveWindowAtts.outputDirectory = "."
  SaveWindowAtts.fileName = filename
  SaveWindowAtts.family = 0
  SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
  SaveWindowAtts.width = 1024
  SaveWindowAtts.height = 1024
  SaveWindowAtts.screenCapture = 0
  SaveWindowAtts.saveTiled = 0
  SaveWindowAtts.quality = 80
  SaveWindowAtts.progressive = 0
  SaveWindowAtts.binary = 0
  SaveWindowAtts.stereo = 0
  SaveWindowAtts.compression = SaveWindowAtts.PackBits  # None, PackBits, Jpeg, Deflate
  SaveWindowAtts.forceMerge = 0
  SaveWindowAtts.resConstraint = SaveWindowAtts.ScreenProportions  # NoConstraint, EqualWidthHeight, ScreenProportions
  SaveWindowAtts.advancedMultiWindowSave = 0
  visit.SetSaveWindowAttributes(SaveWindowAtts)
  visit.SaveWindow()


# either load parameters from file, or enter them afresh using serialize/pickle
# LocalSettings
localVisItDir = "/home/jonny/visit2_10_0.linux-x86_64"
localPythonPackageDir = "/home/jonny/visit2_10_0.linux-x86_64/2.10.0/linux-x86_64/lib/site-packages" 
runRemotely=1
remoteTimeSeriesAstraDB='phase2.wonder.hartree.stfc.ac.uk:/gpfs/stfc/local/HCP084/bwm06/shared/testastra/working/test_analyzeBeam.h5'
remoteTimeSeriesEleSigmaDB='phase2.wonder.hartree.stfc.ac.uk:/gpfs/stfc/local/HCP084/bwm06/shared/testele/phase100/clara.sig-page1.h5'
remoteElePhaseSpaceDB='phase2.wonder.hartree.stfc.ac.uk:/gpfs/stfc/local/HCP084/bwm06/shared/testele/clara_track_electrons_*.vsh5 database'
#remoteElePhaseSpaceDB='phase2.wonder.hartree.stfc.ac.uk:/gpfs/stfc/local/HCP084/bwm06/shared/testele/clara_track_electrons_1.vsh5'
#Remote settings
def isOpenDatabase(filename):
  pass
  return 0

def PhaseSpace(x1,x2,dumpNo):
   visit.AddWindow()
   ScatAttrs = visit.ScatterAttributes()
   ScatAttrs.var1 = x1
   ScatAttrs.var2 = x2
#   ScatAttrs.var1Role = visit.ScatterAtts.Coordinate0
#   ScatAttrs.var2Role = visit.ScatterAtts.Coordinate1
#   ScatAttrs.var3Role = visit.ScatterAtts.None
#   ScatAttsr.var4Role = visit.ScatterAtts.None
   ScatAttrs.SetVar1Role(0) #coord 0
   ScatAttrs.SetVar2Role(1) #coord 1
   ScatAttrs.SetVar3Role(4) #none
   ScatAttrs.SetVar4Role(4) #none   
   ScatAttrs.scaleCube = 0
   visit.AddPlot('Scatter',x1,1,1)
   visit.SetPlotOptions(ScatAttrs)
   visit.SetTimeSliderState(dumpNo)
   visit.DrawPlots() 

def PhaseSpaceXY(elementNo=1):
  if isOpenDatabase(''):
    PhaseSpace('electrons_x','electrons_y',1)
  else:
#    psData=visit.OpenDatabase(remoteElePhaseSpaceDB,elementNo,'Vs')
    psData=visit.OpenDatabase(remoteElePhaseSpaceDB,0,'Vs')
    PhaseSpace('electrons_x','electrons_y',elementNo)

def TimeSeriesS1():
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
  visit.DrawPlots()

    
import sys
sys.path.insert(0,localPythonPackageDir)
import visit
visit.Launch(vdir=localVisItDir)
if runRemotely:
  p2=SetupPhase2()
#  data2=visit.OpenDatabase(remoteTimeSeriesAstraDB,0,'Vs')
  data2=visit.OpenDatabase(remoteTimeSeriesEleSigmaDB,0,'Vs')

else:
  data=visit.OpenDatabase
TimeSeriesS1()
PhaseSpaceXY(2)
SaveWindow('phase2.png')
PhaseSpaceXY(3)
SaveWindow('phase3.png')
PhaseSpaceXY(4)
SaveWindow('phase4.png')

visit.OpenCLI()
visit.OpenGUI()
