#pragma once
// generated during build, do not edit.
#include "include/colour/aces.h"
#include "include/colour/adobergb.h"
#include "include/colour/ergb.h"
#include "include/colour/rec709.h"
#include "include/colour/srgb.h"
#include "include/colour/xyz.h"
#define colour_input_to_xyz colour_ergb_to_xyz
#define colour_xyz_to_input colour_xyz_to_ergb
#define colour_input_print_info colour_ergb_print_info
#define colour_output_to_xyz colour_adobergb_to_xyz
#define colour_xyz_to_output colour_xyz_to_adobergb
#define colour_output_print_info colour_adobergb_print_info
#define colour_camera_to_xyz colour_xyz_to_xyz
#define colour_xyz_to_camera colour_xyz_to_xyz
#define colour_camera_print_info colour_xyz_print_info
