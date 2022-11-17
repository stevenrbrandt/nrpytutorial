import indexedexp as ixp
import reference_metric as rfm
import outputC as outC

# Find the appropriate timestep for the CFL condition.
def add_to_Cfunction_dict_find_timestep():
    # Compute proper distance in all 3 directions.
    delxx = ixp.declarerank1("dxx", DIM=3)
    ds_drn = rfm.ds_dirn(delxx)

    ds_dirn_h = outC.outputC([ds_drn[0], ds_drn[1], ds_drn[2]], ["ds_dirn0", "ds_dirn1", "ds_dirn2"],"returnstring")

    desc="Find the CFL-constrained timestep"
    outC.add_to_Cfunction_dict(
        desc     =desc,
        c_type   ="REAL",
        name     ="find_timestep",
        params   ="const paramstruct *restrict params, REAL *restrict xx[3]",
        preloop  ="REAL dsmin = 1e38; // Start with a crazy high value... close to the largest number in single precision.",
        body     ="REAL ds_dirn0, ds_dirn1, ds_dirn2;\n"+ds_dirn_h+"""
#ifndef MIN
#define MIN(A, B) ( ((A) < (B)) ? (A) : (B) )
#endif
        // Set dsmin = MIN(dsmin, ds_dirn0, ds_dirn1, ds_dirn2);
        dsmin = MIN(dsmin,MIN(ds_dirn0,MIN(ds_dirn1,ds_dirn2)));
""",
        loopopts ="InteriorPoints,Read_xxs,DisableOpenMP",
        postloop ="return dsmin*CFL_FACTOR/wavespeed;\n")


def out_timestep_func_to_file(outfile):
    add_to_Cfunction_dict_find_timestep()
    with open(outfile, "w") as file:
        file.write(outC.outC_function_dict["find_timestep"])

