/********************************************************************/
/*
   Create a dataset with a string in it.
*/
/********************************************************************/

#include <stdlib.h>
#include <string.h>
#include "hdf5.h"
#define FILE "clara-cen.h5"

main() {

  hid_t       file_id, dataset_id, dataspace_id;
  hid_t       new_dataspace_id, new_dataset_id;  /* identifiers */
  hid_t       new_att_id, new_aspace_id;
  herr_t      status;
  char        ** rdata;
  hid_t       dtype, atype;
  size_t      size;
  int         ndims, dim_status, ndimptr, i, maxlen=0;
//  file_id = H5Fcreate(FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  file_id = H5Fopen(FILE, H5F_ACC_RDWR, H5P_DEFAULT);
  dataset_id = H5Dopen(file_id, "/page1/columns/ElementName");
  dataspace_id = H5Dget_space(dataset_id); 
  //  herr_t H5Sselect_all( hid_t dspace_id )
//  dataspace_id = H5Screate (H5S_SCALAR);
  ndims=H5Sget_simple_extent_ndims(dataspace_id);
  printf("ndims: %i\n",ndims);
  hsize_t dimsize[ndims];
  hsize_t maxdims[ndims];
  status=H5Sget_simple_extent_dims(dataspace_id, dimsize,maxdims);
  for (i=0;i<ndims;i++) {
    printf("woohoo\n");
    printf("H5Sget_simple_extentdims: %i\n",dimsize[i]);
  }
  rdata = (char **) malloc (dimsize[0] * sizeof (char *));
  //H5Sget_simple_extent_ndims

  dtype = H5Tcopy (H5T_C_S1);
  atype = H5Tcopy (H5T_C_S1);
  //adata = ""
  // should test we have variable length data in here
  status = H5Tset_size (dtype, H5T_VARIABLE);
  status = H5Dread(dataset_id,dtype,H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata);
  size = 32;
  dtype = H5Tcopy (H5T_C_S1);
  status = H5Tset_size (dtype, size);
  printf("status of H5Tset_size: %i",status);
  status = H5Tset_cset(dtype, H5T_CSET_ASCII);
  new_dataspace_id = H5Screate_simple( ndims, dimsize, maxdims );
  //  space = H5Screate_simple (1, dimsize, NULL);
  printf ("H5Screate_simple returns: %i\n", dataspace_id);
   
  new_dataset_id = H5Dcreate(file_id, "/page1/columns/ElementNameFixedLength", dtype, new_dataspace_id, H5P_DEFAULT);
  status = H5Tset_size (atype,0);
  status = H5Tset_cset (atype, H5T_CSET_ASCII);
  new_aspace_id=H5Screate(H5S_SCALAR);
  new_att_id=H5Acreate(new_dataset_id,"units",atype,new_aspace_id,H5P_DEFAULT);
  char wdata[dimsize[0]][size];
  if (ndims==1) {
    printf("1D data\n");
    for (i=0;i<dimsize[0];i++) {
      printf("allocating\n");
      printf("allocated\n");
      if (strlen(rdata[i])>maxlen) {
        maxlen=strlen(rdata[i]);
      }
      printf("%i: %s, (length %i) \n",i,rdata[i],maxlen);
      printf("strlen %i",strlen(rdata[i]));
      strcpy(wdata[i],rdata[i]);
      printf("in wdata: %s\n",wdata[i]);       
      printf("strlen wdata: %i",strlen(wdata[i]));
    }
    printf("H5WRITING...\n");
    status = H5Dwrite (new_dataset_id, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata);
    printf("status: %i",status);
  }
  else {
    printf("ndims = %i is not supported");
  }
  status = H5Dclose(dataset_id);
  status = H5Sclose(dataspace_id);
  status = H5Fclose(file_id);
}
