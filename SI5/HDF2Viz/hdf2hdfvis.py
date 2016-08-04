# The script has to be run in same directory as your HDF files are.
import numpy
import tables
import re
import sys
import os

if len(sys.argv)==2:
   file_name_in=sys.argv[1]
   print 'Processing file:', file_name_in
else:
   print 'Usage: hdf2hdfvis <FileName> \n'
   sys.exit(1)  
file_name_base  = (file_name_in.split('.')[0]).strip()



input_file=tables.open_file(file_name_in,'r')






# Below lines dump HDF5 strucutre into files and then change its shape for accessing tables
file_output=open('H5LS_OUT.TEXT','w')
sys.stdout=file_output
print(input_file)
sys.stdout = sys.__stdout__
file_output.close()


# Convert dumped HDF5 fil   e structure to necessary arrays call
help_array=[]
hand = open('H5LS_OUT.TEXT')
# print hand
#wyjsciowka=[]
help_file=open('HDF5STRUCT.TEXT','w')
add_path_data_beginning='input_file.root'
add_read_data_end='.read()'
for line in hand:
	line = line.rstrip()
	x=re.findall('\S+columns\S+',line)
# This loop checks reomoves or changes unnecessary strings in file which contatins structure of HDF source file.
	if len(x) > 0 :
            sx=str(x)
            sx=sx.replace("[","")
            sx=sx.replace("]","")
            sx=sx.replace("'","")
            sx=sx.replace("/",".")
            sx=add_path_data_beginning+sx+add_read_data_end
            help_array.append(sx)
            # print sx



help_file.write("\n".join(help_array))

help_file.close()
# End of array structure processing, file left for further diagnostics and purpose


# Load data from source file basing on above detected file structure (directories).

#input_file=tables.open_file(file_name_in,'r')

#Read names from file generated by h5ls -rf filename.hdf > output.txt
read_names=open('HDF5STRUCT.TEXT','r')
# Read the number of columns and set counter of columns
counter=len(open('HDF5STRUCT.TEXT','r').readlines())
# print counter
# Read of names of datasets in columns
array_items=read_names.read().splitlines()



output_file=tables.open_file(file_name_base+'VIS.hdf5','w')
working_array=numpy.array([eval(array_items[0])])



for counter2 in range (1,counter):
    # print array_items[counter2]
    working_array=numpy.append(working_array,numpy.array([eval(array_items[counter2])]),axis=0)


working_array=working_array.transpose()

# print working_array
ParticleGroup=output_file.create_array('/','Particles', working_array)



# print array_items
runGroup=output_file.create_group('/','runInfo','')
runGroup._v_attrs.vsType='runInfo'
boundsGroup=output_file.create_group('/','globalGridGlobalLimits','')
boundsGroup._v_attrs.vsType='limits'
boundsGroup._v_attrs.vsKind='Cartesian'
# This generates lower and upper bounds for VizSchema but is not necessary (?)
#boundsGroup._v_attrs.vsLowerBounds=numpy.array([numpy.min(column_x),numpy.min(column_y)])
#boundsGroup._v_attrs.vsUpperBounds=numpy.array([numpy.max(column_x),numpy.max(column_y)])

# Passes names of all columns to VizSchema HDF5

array_items[0]=re.sub('\S+.columns.','',array_items[0])
array_items[0]=re.sub('.read\S+','',array_items[0])

for counter2 in range (1,counter):
    array_items[counter2]=re.sub('\S+.columns.',',',array_items[counter2])
    array_items[counter2]=re.sub('.read\S+','',array_items[counter2])
    array_items[counter2]=str(array_items[counter2])
ParticleGroup._v_attrs.vsLabels=''.join(array_items)

# VizSchema required strings - need to discuss it with Jonathan.
timeGroup=output_file.create_group('/','time','time')
timeGroup._v_attrs.vsType='time'
ParticleGroup._v_attrs.vsType='variableWithMesh'
ParticleGroup._v_attrs.vsTimeGroup='time'
# Checks the number of dimensions and then put proper number in file.
if counter==1:
    ParticleGroup._v_attrs.vsNumSpatialDims = 1
elif counter==2:
    ParticleGroup._v_attrs.vsNumSpatialDims = 2
elif counter>=3:
    ParticleGroup._v_attrs.vsNumSpatialDims = 3
else:
    ParticleGroup._v_attrs.vsNumSpatialDims = 0
ParticleGroup._v_attrs.vsLimits='globalGridGlobalLimits'
os.remove("H5LS_OUT.TEXT")
os.remove("HDF5STRUCT.TEXT")
output_file.close()