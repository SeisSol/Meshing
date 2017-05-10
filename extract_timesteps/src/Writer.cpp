#include "Writer.h"
#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>

static void check_err(int const stat, int const line, char const* file)
{
  if (stat != NC_NOERR) {
    fprintf(stderr,"line %d of %s: %s\n", line, file, nc_strerror(stat));
    fflush(stderr);
    exit(1);
  }
}

void writeTimesteps(std::vector<double> const& timesteps, std::string const& fileName)
{
  int  stat;  /* return status */
  int  ncid;  /* netCDF id */

  /* dimension ids */
  int element_dim;

  /* dimension lengths */
  size_t element_len = timesteps.size();

  /* variable ids */
  int timestep_id;

  /* rank (number of dimensions) for each variable */
#   define RANK_timestep 1

  /* variable shapes */
  int timestep_dims[RANK_timestep];

  /* enter define mode */
  stat = nc_create(fileName.c_str(), NC_CLOBBER, &ncid);
  check_err(stat,__LINE__,__FILE__);

  /* define dimensions */
  stat = nc_def_dim(ncid, "element", element_len, &element_dim);
  check_err(stat,__LINE__,__FILE__);

  /* define variables */

  timestep_dims[0] = element_dim;
  stat = nc_def_var(ncid, "timestep", NC_DOUBLE, RANK_timestep, timestep_dims, &timestep_id);
  check_err(stat,__LINE__,__FILE__);

  /* assign per-variable attributes */
  stat = nc_put_att_text(ncid, timestep_id, "units", 1, "s");
  check_err(stat,__LINE__,__FILE__);


  /* leave define mode */
  stat = nc_enddef (ncid);
  check_err(stat,__LINE__,__FILE__);

  /* assign variable data */
  stat = nc_put_var_double(ncid, timestep_id, &timesteps.data()[0]);
  check_err(stat,__LINE__,__FILE__);

  stat = nc_close(ncid);
  check_err(stat,__LINE__,__FILE__);
}
