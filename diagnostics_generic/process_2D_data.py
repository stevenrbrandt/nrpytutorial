#
# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com
import numpy as np
from scipy.interpolate import griddata


def extract_1D_slice_from_2D_data(datafile, y_output_value,
                                  col_number_x, col_number_y, col_number_data,
                                  xminmax, sample_numpts_x=100, interp_method='linear'):
    cols_list = np.loadtxt(datafile).T # Transposed for easier unpacking

    input_x = cols_list[col_number_x]
    input_y = cols_list[col_number_y]
    input_data = cols_list[col_number_data]

    num_pts = len(input_x)

    # Define input grid
    input_points = np.zeros((num_pts, 2))
    for i in range(num_pts):
        input_points[i][0] = input_x[i]
        input_points[i][1] = input_y[i]

    # Define output grid
    grid_size_x = complex(0, sample_numpts_x)
    grid_size_y = complex(0, 1)
    pl_xmin = xminmax[0]
    pl_xmax = xminmax[1]
    pl_ymin = y_output_value
    pl_ymax = y_output_value
    output_grid_x, output_grid_y = np.mgrid[pl_xmin:pl_xmax:grid_size_x, pl_ymin:pl_ymax:grid_size_y]

    # Interpolate input grid data to output grid
    output_grid_data = griddata(input_points, input_data, (output_grid_x, output_grid_y), method=interp_method)
    return output_grid_x, output_grid_data


def generate_uniform_2D_grid(datafile,
                             col_number_x, col_number_y, col_number_data,
                             xminmax, yminmax,
                             sample_numpts_xy=100, interp_method='cubic'):
    cols_list = np.loadtxt(datafile).T  # Transposed for easier unpacking

    input_x = cols_list[col_number_x]
    input_y = cols_list[col_number_y]
    input_data = cols_list[col_number_data]

    num_pts = len(input_x)

    # Define input grid
    input_points = np.zeros((num_pts, 2))
    for i in range(num_pts):
        input_points[i][0] = input_x[i]
        input_points[i][1] = input_y[i]

    # Define output grid
    grid_size_x = complex(0, sample_numpts_xy)
    grid_size_y = complex(0, sample_numpts_xy)
    pl_xmin = xminmax[0]
    pl_xmax = xminmax[1]
    pl_ymin = yminmax[0]
    pl_ymax = yminmax[1]
    output_grid_x, output_grid_y = np.mgrid[pl_xmin:pl_xmax:grid_size_x, pl_ymin:pl_ymax:grid_size_y]

    # Interpolate input grid data to output grid
    output_grid_data = griddata(input_points, input_data, (output_grid_x, output_grid_y), method=interp_method)

    return output_grid_x, output_grid_y, output_grid_data.T  # Transpose back to original (see above)

