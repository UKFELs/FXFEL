import os, sys, glob
import tables, numpy
import optparse
import subprocess
######
# (c) Jonathan Smith, Tech-X UK Ltd 2015
######
"""
This script loops through a directory, and takes any sdds files it can find and runs sdds2hdf on them
With the output, it tests whether this is a particle data set of the type generated by Elegant
and if so adds appropriate metadata to allow visualisation in VisIt 
"""

datadir="."

def filenameToBasename(filename):
  basename=filename[:-3]
  return basename

def isSDDSelegantParticleData(someData): 
  columns=someData.columns._v_children.keys()
  print columns
  if ('x' in columns and 'y' in columns and 't' in columns and 'xp' in columns and 'yp' in columns and 'dt' in columns): 
    if 'Step' in someData.parameters._v_children.keys():
      if 's' in someData.parameters._v_children.keys():
        if 'medt' in someData.parameters._v_children.keys():
          print 'median time: '+str(someData.parameters.medt[0])+str(someData.parameters.medt._v_attrs.units)
        print 'time from s: '+str(float(someData.parameters.s[0]/2.998e8))
        print 's: '+str(someData.parameters.s[0])+str(someData.parameters.s._v_attrs.units)
        return 1
      else:
        print "WARNING: No s value"
        return 0
    else:
      print "WARNING: No Step Information"
      return 0
  else:
    print "WARNING: Not the correct columns for particle phase space"
    return 0

def sddsh5ptclstoVsH5(filename):
  h5=tables.open_file(filename)
# Columns on top level or pages - always pages I think. This needs expansion
# to handle the multiple pages case, as happens in parameter sweeps, etc.
  for page in h5.root._f_list_nodes():
    print page._v_name
    sddsH5Data=page
    if 'columns' in sddsH5Data._v_children.keys():
  #  if sddsH5Data._v_children:
        if isSDDSelegantParticleData(sddsH5Data):
  # Write a vizSchema file
          params=sddsH5Data.parameters
          for i in params._v_children.keys():
            print i+': '+str(params._f_get_child(i).read())+' '+str(params._f_get_child(i)._v_attrs.units)          
          writeSDDSptclsToVsH5(sddsH5Data,filename[:-3],'electrons')
        else:
          print sddsH5Data._v_children    
          print "Not a particle data set"
  h5.close()

    
def writeSDDSptclsToVsH5(sddsH5Data, basename, particleName):
    vsh5=tables.open_file(basename+'_'+particleName+'.vsh5','w')
    if 'particleID' in sddsH5Data.columns._v_children.keys():
      electrons=numpy.array([sddsH5Data.columns.x.read(), sddsH5Data.columns.y.read(), sddsH5Data.columns.t.read(), sddsH5Data.columns.xp.read(), sddsH5Data.columns.yp.read(), sddsH5Data.columns.dt.read(), sddsH5Data.columns.p.read(), sddsH5Data.columns.particleID.read()]).transpose()
    else:
      electrons=numpy.array([sddsH5Data.columns.x.read(), sddsH5Data.columns.y.read(), sddsH5Data.columns.t.read(), sddsH5Data.columns.xp.read(), sddsH5Data.columns.yp.read(), sddsH5Data.columns.dt.read(), sddsH5Data.columns.p.read()]).transpose()     
    print electrons
    elecGroup=vsh5.create_array('/','electrons',electrons)
    runGroup=vsh5.create_group('/','runInfo','')
    runGroup._v_attrs.vsType='runInfo'
    boundsGroup=vsh5.create_group('/','globalGridGlobalLimits','')
    boundsGroup._v_attrs.vsType='limits'
    boundsGroup._v_attrs.vsKind='Cartesian'
    boundsGroup._v_attrs.vsLowerBounds=numpy.array([numpy.min(sddsH5Data.columns.x.read()),numpy.min(sddsH5Data.columns.y.read()),numpy.min(sddsH5Data.columns.t.read())])
    boundsGroup._v_attrs.vsUpperBounds=numpy.array([numpy.max(sddsH5Data.columns.x.read()),numpy.max(sddsH5Data.columns.y.read()),numpy.max(sddsH5Data.columns.t.read())])
    timeGroup=vsh5.create_group('/','time','actually, position s')
    timeGroup._v_attrs.vsStep=sddsH5Data.parameters.Step.read()
    timeGroup._v_attrs.vsType='time'
    elecGroup._v_attrs.vsType='variableWithMesh'
    elecGroup._v_attrs.vsTimeGroup='time'
    elecGroup._v_attrs.vsNumSpatialDims = 3
    elecGroup._v_attrs.vsLimits='globalGridGlobalLimits'
    elecGroup._v_attrs.vsLabels=particleName+'_x,'+particleName+'_y,'+particleName+'_t,'+particleName+'_xp,'+particleName+'_yp,'+particleName+'_dt,'+particleName+'_p,particleID'
    elecGroup._v_attrs.vsUnits='m,m,s,,,,$m_{e}c$,'
    varsGroup=vsh5.create_group('/',particleName+'_z','')
    varsGroup._v_attrs.vsType='vsVars'
    varsGroup._v_attrs.z='2.998e8*'+particleName+'_t'
    varsGroup=vsh5.create_group('/',particleName+'_tmod','')
    varsGroup._v_attrs.vsType='vsVars'
    varsGroup._v_attrs.tmod='('+particleName+'_t-'+str(numpy.average(sddsH5Data.columns.t.read()))+")"
    varsGroup=vsh5.create_group('/','zmod','')
    varsGroup._v_attrs.vsType='vsVars'
    varsGroup._v_attrs.zmod='2.998e8*'+particleName+'_tmod'
    varsGroup=vsh5.create_group('/','electrons_dtmod','')
    varsGroup._v_attrs.vsType='vsVars'
    varsGroup._v_attrs.zmod=particleName+'_dt-'++str(numpy.average(sddsH5Data.columns.dt.read()))
    if 's' in sddsH5Data.parameters._v_children.keys():
      elecGroup._v_attrs.time=sddsH5Data.parameters.s.read()
      timeGroup._v_attrs.vsTime=sddsH5Data.parameters.s.read()
    for paramName in sddsH5Data.parameters._v_children.keys():
      runGroup._f_setattr('vs'+paramName,str(sddsH5Data.parameters._f_get_child(paramName)[0])+sddsH5Data.parameters._f_get_child(paramName)._v_attrs.units)
    vsh5.close()
#os.chdir(r'c:\Users\Jonny\Documents\SDDS')


#os.chdir(datadir)
#sddsFiles = glob.glob('*.sdds')
#sddsCenFiles = glob.glob('*.cen')
#sddsTwiFiles = glob.glob('*.twi')
#sddsSigFiles = glob.glob('*.sig')
#sddsH5Files = glob.glob('*.h5')
#print sddsFiles
#for sddsFile in sddsFiles:
if (len(sys.argv) == 2):
#  sddsFile
#  print sddsFile
#  h5fileName=sddsFile[:-4]+'h5'
  h5fileName=sys.argv[1]
  print "using "+h5fileName+" as input file"
  if os.path.isfile(os.getcwd()+os.sep+h5fileName):
    sddsh5ptclstoVsH5(h5fileName)
    try:
      pass
    except:
      print "Failed to perform conversion"
  else:
    print "Didn't find h5 file converted directly from sdds"
    try:
      print "Writing h5 file using sdds2hdf"
      stdoutput=os.system(sdds2hdfCommand+' '+sddsFile+' '+h5fileName)
      print stdoutput
      try:
        sddsh5ptclstoVsH5(h5fileName)
      except:
        print "Failed to perform conversion"
      
    except:
      print "Failed to write h5 from sdds2hdf"
# r'C:\Users\Jonny\Downloads\01_250-MC-08.0350.001\01_250-MC-08.0350.001\CLA-MU1-DIA-SCR-01.h5'  
# Something wrong with subprocess in my python distribution... 
# below would be more contemporary as os.system is getting phased out
#  stdoutput=subprocess.check_output([args])
#print sddsH5Files
