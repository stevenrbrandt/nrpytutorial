#
# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

# Output data on the plane closest to xy or xz
def output_plane_xz_or_xy_body(plane="xz", include_ghosts=True):
    out_str = r"""  // First unpack:
  const paramstruct params = griddata->params;
  REAL *restrict xx[3] = { griddata->xx[0],griddata->xx[1],griddata->xx[2] };
  const int Nx0 = params.Nxx0, Nx1 = params.Nxx1, Nx2 = params.Nxx2;
  const int Nxx_plus_2NGHOSTS0 = params.Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params.Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params.Nxx_plus_2NGHOSTS2;

  // Next set loop bounds, based on coordinate system (params.CoordSystemName)
"""
    if include_ghosts:
        out_str += r"""
  int numpts_i0=Nx0 + 2*NGHOSTS, numpts_i1=Nx1 + 2*NGHOSTS, numpts_i2=Nx2 + 2*NGHOSTS;

  int i0_pts[Nx0 + 2*NGHOSTS], i1_pts[Nx1 + 2*NGHOSTS], i2_pts[Nx2 + 2*NGHOSTS];
#pragma omp parallel for
  for(int i0=0; i0<Nx0 + 2*NGHOSTS; i0++) i0_pts[i0] = i0;
#pragma omp parallel for
  for(int i1=0; i1<Nx1 + 2*NGHOSTS; i1++) i1_pts[i1] = i1;
#pragma omp parallel for
  for(int i2=0; i2<Nx2 + 2*NGHOSTS; i2++) i2_pts[i2] = i2;
"""
    else:
        out_str += r"""
  int numpts_i0=Nx0, numpts_i1=Nx1, numpts_i2=Nx2;

  int i0_pts[Nx0], i1_pts[Nx1], i2_pts[Nx2];
#pragma omp parallel for
  for(int i0=NGHOSTS; i0<Nx0 + NGHOSTS; i0++) i0_pts[i0-NGHOSTS] = i0;
#pragma omp parallel for
  for(int i1=NGHOSTS; i1<Nx1 + NGHOSTS; i1++) i1_pts[i1-NGHOSTS] = i1;
#pragma omp parallel for
  for(int i2=NGHOSTS; i2<Nx2 + NGHOSTS; i2++) i2_pts[i2-NGHOSTS] = i2;
"""
    if plane == "xz":
        out_str += r"""
  if(strstr(params.CoordSystemName, "Cartesian") != NULL) {
    // xz-plane == { y_mid }, where y index is i1
    numpts_i1 = 1;
    i1_pts[0] = (Nx1 + 2*NGHOSTS) / 2;
  } else if(strstr(params.CoordSystemName, "Cylindrical") != NULL) {
    // xz-plane == { phi_min or phi_mid }, where phi index is i1; note that phi_min=-PI and phi_max=+PI, modulo ghostzones
    numpts_i1 = 2;
    i1_pts[0] = NGHOSTS;
    i1_pts[1] = (Nx1 + 2*NGHOSTS) / 2;
  } else if(strstr(params.CoordSystemName, "Spherical") != NULL ||
            strstr(params.CoordSystemName, "SymTP")     != NULL) {
    // xz-plane == { phi_min or phi_mid }, where phi index is i2; note that phi_min=-PI and phi_max=+PI, modulo ghostzones
    numpts_i2 = 2;
    i2_pts[0] = NGHOSTS;
    i2_pts[1] = (Nx2 + 2*NGHOSTS) / 2;
  } else {
    fprintf(stderr, "output_"""+plane+r"""...() ERROR: params.CoordSystemName == %s unsupported.\n", params.CoordSystemName);
    exit(1);
  }
"""
    elif plane == "xy":
        out_str += r"""
  if(strstr(params.CoordSystemName, "Cartesian") != NULL) {
    // xy-plane == { z_mid }, where z index is i2
    numpts_i2 = 1;
    i2_pts[0] = (Nx2 + 2*NGHOSTS) / 2;
  } else if(strstr(params.CoordSystemName, "Cylindrical") != NULL) {
    // xy-plane == { z_mid }, where z index is i2
    numpts_i2 = 1;
    i2_pts[0] = (Nx2 + 2*NGHOSTS) / 2;
  } else if(strstr(params.CoordSystemName, "Spherical") != NULL ||
            strstr(params.CoordSystemName, "SymTP")     != NULL) {
    // xy-plane == { theta_mid }, where theta index is i1
    numpts_i1 = 1;
    i1_pts[0] = (Nx1 + 2*NGHOSTS) / 2;;
  } else {
    fprintf(stderr, "output_"""+plane+r"""...() ERROR: params.CoordSystemName == %s unsupported.\n", params.CoordSystemName);
    exit(1);
  }
"""
    return out_str
