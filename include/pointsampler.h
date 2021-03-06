/*
    This file is part of corona-13.

    corona-13 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    corona-13 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with corona-13. If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include "mf.h"
#include <stdint.h>
#include <stdio.h>


struct pointsampler_t;
typedef struct pointsampler_t pointsampler_t;

// fwd declare
struct path_t;

void pointsampler_print_info(FILE *f);

pointsampler_t *pointsampler_init(uint64_t frame);

// get random number:
float pointsampler(struct path_t *p, int dim);

void pointsampler_splat(struct path_t *p, mf_t value);

// return 1 if the new path (tentative) is accepted over the current one.  will
// also call view_splat() to accumulate into the frame buffer (only upon
// rejection)
int pointsampler_accept(struct path_t *curr, struct path_t *tent);

//Set the factor for the spp for the designated pixelblock
void setBlockSamples(int blockNumber, double value);
//disable welch sampling and enable factored sampling according to the variance estimation
void enableFactoredSampling();
//Return the factor of the pixel block of the given pixels (gets called in blackmanharris.h)
int getFactor(float i, float j);
//The given samples are returned as a .pfm file
void write_samples_as_framebuffer();
//Return the sum of all samples that are distributed during the next iteration
uint64_t getSampleCount(int width, int height);

// mutate path into new_path
void pointsampler_mutate(struct path_t *curr, struct path_t *tent);
void pointsampler_mutate_with_pixel(struct path_t *curr, struct path_t *tent, float i, float j);

void pointsampler_clear();
void pointsampler_cleanup(pointsampler_t *s);
void pointsampler_set_large_step(pointsampler_t *t, float p_large_step);

// accum rest of current sample, in case it was never rejected before the end
// of the render.
void pointsampler_finalize(pointsampler_t *s);
// called before every progression
void pointsampler_prepare_frame(pointsampler_t *s);

// interface to inject fake random numbers from behind:
void pointsampler_enable_fake_random(pointsampler_t *s);
void pointsampler_disable_fake_random(pointsampler_t *s);
void pointsampler_set_fake_random(pointsampler_t *s, int dim, float rand);

void pointsampler_stop_learning(pointsampler_t* s);

