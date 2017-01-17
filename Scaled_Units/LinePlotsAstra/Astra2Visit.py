#!/usr/bin/python
import glob,sys,os,re,tables,numpy
if len(sys.argv) < 2:
  print "Usage : Astra2Visit basename"
  exit()
basename=sys.argv[1]
print "basename "+basename
allfilelist=glob.glob(basename+'*')
print allfilelist
particleFileList=[]
positions={}
for fname in allfilelist:
  mg=re.match(basename+'\.\d+\.(\d+)\.\d+',fname)
  if mg:
    pos=float(mg.groups()[0])/100.
    particleFileList.append(fname)
    print pos
    positions[fname]=pos

print "particle files\n\n"
print particleFileList

MachinePosS=[]
DataSets={}
DataSetVals={}
for fname in sorted(positions, key=positions.get):
#particleFileList:
  try: 
    os.system('astra2elegant '+fname+' '+fname+'.sdds')
  except:
    print "conversion to sdds from astra failed"
  try:
    os.system('sddsanalyzebeam '+fname+'.sdds '+fname+'.analyze.sdds')
    os.system('sdds2hdf '+fname+'.analyze.sdds '+fname+'.analyze.h5')  
  except:
    print "something else failed"
    
  mg=re.match(basename+'\.\d+\.(\d+)\.\d+',fname)
  if mg:
    pos=float(mg.groups()[0])/100.
    print pos
  if pos in MachinePosS:
    print "position "+str(pos)+" from "+fname+" is already recognised in array, do we have duplicate data at some location?"
    exit(2)
  else:
    MachinePosS.append(pos)
    print "Working on "+fname+".analyze.h5..."
    h5=tables.open_file(fname+".analyze.h5","r")
    pages=h5.root._v_children
    for page in pages:
      if page[:4]!="page":
        print "The data structure has no top level page. Looks like it's not from sdds, and so not for this script"
        exit(3)
    else:
      thisFileDatasets=[]
#      vsh5out=tables.open_file(infile[:-3]+'-'+page+".h5","w")
      print "page:",page
      if page=="page2":
        print "Please contact the developers to enable this script to handle multiple pages in your sdds files"
        exit(6)
      for i in h5.walk_nodes("/"+page):
        if type(i)==tables.group.RootGroup:
          print "root: ",i._v_name
          print type(i)
        elif type(i)==tables.unimplemented.UnImplemented:
          print "skipping: ",i._v_name
          print "Please notify the developers that some string has not been converted prior to running this script"
#          exit()
        if type(i)==tables.group.Group:
          print "Group: "+i._v_name
        elif type(i)==tables.array.Array:
          #print i._v_name
          thisFileDatasets.append(i)
          if i._v_name in DataSets:
            #print "Already the dataset exists"
            if i.read().shape==(1,):
              DataSets[i._v_name][str(pos)]=i.read()[0]
              DataSetVals[i._v_name].append(i.read()[0])
            else:
              print "but apparently wrong shape"
              exit(4)
          else:
            if i.read().shape==(1,):
              print "adding "+i._v_name+' dataset, but not checking type'
              #DataSets[i._v_name]={ }
              DataSetVals[i._v_name]=[ i.read()[0] ]
              DataSets[i._v_name]={"name":i._v_name,"type":type(i.read()[0]),str(pos):i.read()[0]}
            else:
              print "but apparently wrong shape"
              exit(5)
        else:
          print type(i)," is not handled"
    h5.close()
    os.remove(fname+'.sdds')
    os.remove(fname+'.analyze.sdds')
    os.remove(fname+'.analyze.h5')

#  print DataSets
print DataSetVals

vsh5out=tables.open_file(basename+"_analyzeBeam.h5","w")
print "Writing Mesh"
vsh5out.create_array('/','spatialPositions', numpy.array(MachinePosS))
vsh5out.root.spatialPositions._v_attrs.vsType="mesh"
vsh5out.root.spatialPositions._v_attrs.vsKind="structured"

for dataset in DataSets:
  print dataset
  print "Recording "+dataset
  vsh5out.create_array('/',dataset,numpy.array(DataSetVals[dataset]))
  vsh5out.root._v_children[dataset]._v_attrs.vsMesh="/spatialPositions"
  vsh5out.root._v_children[dataset]._v_attrs.vsType = "variable"
  #vsh5out.root._v_children[dataset._v_name]._v_attrs.vsUnits=dataset._v_attrs['units']
  print vsh5out.root._v_children[dataset]
  print vsh5out.root._v_children[dataset]._v_attrs.vsMesh

#      vsh5out.
vsh5out.create_group('/','runInfo','')
vsh5out.close()
