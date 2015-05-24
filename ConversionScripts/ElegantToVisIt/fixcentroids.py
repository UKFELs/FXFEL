#!/bin/python
import tables,numpy,sys
infile=sys.argv[1]
print infile
h5=tables.open_file(infile,"r")
pages=h5.root._v_children
for page in pages:
 if page[:4]!="page":
  print "The data structure has no top level page. Looks like it's not from sdds, and so not for this script"
 else:
  datasets=[]
  vsh5out=tables.open_file(infile[:-3]+'-'+page+".h5","w")
  print "page:",page
  for i in h5.walk_nodes("/"+page):
    if type(i)==tables.group.RootGroup:
      print "root: ",i._v_name
      print type(i)
    elif type(i)==tables.unimplemented.UnImplemented:
      print "skipping: ",i._v_name
    if type(i)==tables.group.Group:
      print "Group: "+i._v_name
    elif type(i)==tables.array.Array:
      print i._v_name
      datasets.append(i)
    else:
      print type(i)," is not handled"
#  print datasets
  for dataset in datasets:
    if dataset._v_name=="s":
      print "This is our mesh"     
      vsh5out.create_array('/','spatialPositions', h5.root._v_groups[page].columns.s.read())
      vsh5out.root.spatialPositions._v_attrs.vsType="mesh"
      vsh5out.root.spatialPositions._v_attrs.vsKind="structured"
    else:
      print dataset._v_name
      vsh5out.create_array('/',dataset._v_name,dataset.read())
      vsh5out.root._v_children[dataset._v_name]._v_attrs.vsMesh="/spatialPositions"
      vsh5out.root._v_children[dataset._v_name]._v_attrs.vsType = "variable"
      vsh5out.root._v_children[dataset._v_name]._v_attrs.vsUnits=dataset._v_attrs['units']
      print vsh5out.root._v_children[dataset._v_name]
      print vsh5out.root._v_children[dataset._v_name]._v_attrs.vsMesh

#      vsh5out.
  vsh5out.create_group('/','runInfo','')

  vsh5out.close()
h5.close()
