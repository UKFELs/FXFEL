import sys,os,tables,numpy,re
filename = sys.argv[1]
h5=tables.open_file(filename,'r')
vsh5out=tables.open_file(filename[:-3]+'.vsh5','w')
#make list of unique parameter names
#for each sliceNN of those write a dataset.
datasets=[]
slices=[] 
for param in h5.walk_nodes("/page1/parameters/"):
  pname= param._v_name
  print pname
  if "Slice" in param._v_name:    
    if param._v_name[-7:-2]=="Slice":
      if pname[:-7] in datasets:
        #print "slice :"+pname[:-7]+' : '+pname[-2:]
        print ".",
      else:
        datasets.append(pname[:-7])
      if pname[-2:] in slices:
        print ",",
      else:
        slices.append(pname[-2:])
    else:
      print "Slice numbering is not two characters, so the code around here needs to change as our assumptions are incorrect. Contact a developer"

print datasets
print slices
print numpy.array(slices)
sliceArray=numpy.array(slices).astype(float)
print sliceArray
vsh5out.create_array('/','sliceNo', sliceArray)
vsh5out.root.sliceNo._v_attrs.vsType="mesh"
vsh5out.root.sliceNo._v_attrs.vsKind="structured"
for datasetName in datasets:
  dataset=numpy.zeros(len(slices)) #98 most likely.
  for slice in slices:
# python numbering starts at zero, data was likely creating with fortran numbering of slices starting at 1.
# Also expect each parameter to contain only one value, but stored in an array hence trailing [0]
    dataset[int(slice)-1]=h5.root.page1.parameters._v_children[datasetName+"Slice"+slice][0]
  
  vsh5out.create_array('/',datasetName,dataset)
  vsh5out.root._v_children[datasetName]._v_attrs.vsMesh="/sliceNo"
  vsh5out.root._v_children[datasetName]._v_attrs.vsType = "variable"
#  vsh5out.root._v_children[datasetName]._v_attrs.vsUnits=dataset._v_attrs['units']
vsh5out.close()
h5.close()
