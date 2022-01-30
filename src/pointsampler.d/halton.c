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

#include "corona_common.h"
#include "pointsampler.h"
#include "sampler.h"
#include "render.h"
#include "pathspace.h"
#include "points.h"
#include "threads.h"
#include "ext/halton/halton.h"
#include "framebuffer.h"

#include <stdio.h>
#include <float.h>
#include <pthread.h>

//Number of color channels, standard value is 3
#define NUMCHANNELS 3
//The size of each grid cell for rendering once welch sampling is disabled
#define GRID_SIZE 16
//Maximum amount of samples per pixel
#define MAXSAMPLEVALUE 40
//Samples per Pixel per Pixel Block are saved here
static int *sample_factor;
static int enableFactorSampling = 0;
static int init = 0;
static int gridnumber = 0;
static int spp = 0;
static int row = 0;
static int col = 0;
//10 is a placeholder, value is set in pointsampler_mutate
static int num_horizontal_cells = 10;
static int num_vertical_cells = 10;
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

typedef struct fake_randoms_t
{
  int enabled;
  float rand[40];
}
fake_randoms_t;

typedef struct pointsampler_t
{
  fake_randoms_t *rand;
  halton_t h;
  uint64_t reinit;
}
pointsampler_t;

void pointsampler_print_info(FILE *f)
{
  fprintf(f, "mutations: halton points\n");
}

pointsampler_t *pointsampler_init(uint64_t frame)
{
  pointsampler_t *s = (pointsampler_t *)malloc(sizeof(pointsampler_t));
  s->rand = calloc(rt.num_threads, sizeof(*s->rand));
  s->reinit = 0;
  halton_init_random(&s->h, frame);
  return s;
}

int pointsampler_accept(path_t *curr, path_t *tent) { return 0; }
void pointsampler_clear() {}
void pointsampler_cleanup(pointsampler_t *s)
{
  free(s->rand);
  free(s);
}
void pointsampler_set_large_step(pointsampler_t *t, float p_large_step) {}
void pointsampler_finalize(pointsampler_t *s) {}

float pointsampler(path_t *p, int i)
{
  const int tid = common_get_threadid();
  if(rt.pointsampler->rand[tid].enabled)
    return rt.pointsampler->rand[tid].rand[i];

  int v = p->length;
  const int end = p->v[v].rand_beg;
  const int dim = end + i;
  if(dim >= halton_get_num_dimensions())
    // degenerate to pure random mersenne twister
    return points_rand(rt.points, common_get_threadid());
  else
    // note that this clips the bits in p->index to 32:
    return halton_sample(&rt.pointsampler->h, dim, p->index);
}

void pointsampler_splat(path_t *p, mf_t value)
{
  render_splat(p, value);
}

int getFactor(float i, float j) {
  //If called while welch sampling is still enabled
  if (!enableFactorSampling) return 1;
  int factor = 1;
  //Calculate the factor for the Gridcell, for normalizing
  int row = floor(i / GRID_SIZE);
  int col = floor(j / GRID_SIZE);
  factor = sample_factor[col * num_horizontal_cells + row];
  return factor;
}

void setBlockSamples(int blockNumber, double value) {
  int blockValue = (int)floor(value);
  //clamp to 1 or the max value
  if(blockValue < 1) {
    blockValue = 1;
  } else if(blockValue > MAXSAMPLEVALUE) {
    blockValue = MAXSAMPLEVALUE;
  } else {
    sample_factor[blockNumber] = blockValue;
  }
  //printf("Set factor of Block %d to %d \n", blockNumber, blockValue);
}

void enableFactoredSampling() { 
  enableFactorSampling = 1;
  write_samples_as_framebuffer();
}

void write_samples_as_framebuffer() {
  framebuffer_t *fb = malloc(sizeof(framebuffer_t));
  framebuffer_header_t *fb_header = malloc(sizeof(framebuffer_header_t));

  fb_header->channels = 3;
  fb_header->gain = 1;
  fb_header->height = num_vertical_cells;
  fb_header->width = num_horizontal_cells;

  fb->header = fb_header;
  fb->retain = 0;

  //divide by the highest possible value so the resulting values are between 0 and 1
  float* buffer = (float *) malloc(3 * num_horizontal_cells*num_vertical_cells * sizeof(float));
  for (int i = 0; i < num_horizontal_cells * num_vertical_cells; i++) {
    buffer[3*i] = (sample_factor[i] + 0.0f) / MAXSAMPLEVALUE;
    buffer[3*i+1] = (sample_factor[i] + 0.0f) / MAXSAMPLEVALUE;
    buffer[3*i+2] = (sample_factor[i] + 0.0f) / MAXSAMPLEVALUE;

  }
  fb->fb = buffer;

  char* name = "SampleDistribution.pfm";
  fb_export(fb, name, 0, 0);
  free(fb);
  free(fb_header);
  free(buffer);
}

void pointsampler_mutate(path_t *curr, path_t *tent)
{
  //init global variables
  if (!init) {
    init = 1;
    sample_factor = (int*) malloc(view_width() / GRID_SIZE * view_height() / GRID_SIZE * sizeof(int));
    num_horizontal_cells = view_width() / GRID_SIZE;
    num_vertical_cells = view_height() / GRID_SIZE;
    for (int i = 0; i < num_horizontal_cells*num_vertical_cells; i++) {
      sample_factor[i] = 1;
    }
  }

  if(!enableFactorSampling) {
    //standard procedure
    path_init(tent, tent->index, tent->sensor.camid);
    sampler_create_path(tent);
  } else {
    //calculate random value within given pixel
    float i = points_rand(rt.points, common_get_threadid());
    i = i + row + (gridnumber % num_horizontal_cells) * GRID_SIZE;
    //float j = pointsampler(curr, s_dim_image_y);
    //printf("%f is j \n", j);
    float j = points_rand(rt.points, common_get_threadid());
    j = j + col + (gridnumber / num_horizontal_cells) * GRID_SIZE;
    pointsampler_mutate_with_pixel(curr, tent, i, j);

    pthread_mutex_lock(&mutex);
    //New sample within same pixel
    ++spp;
    //new Pixel (row gets incremented)
    if (spp >= sample_factor[gridnumber] * NUMCHANNELS) {
      spp = 0;
      ++row;
      //New Pixel (col gets incremented)
      if (row >= GRID_SIZE) {
        ++col;
        row = 0;
        //New Grid cell
        if (col >= GRID_SIZE) {
          col = 0;
          ++gridnumber;
          //Change into first grid cell since all grid cells have been sampled
          if (gridnumber >= num_horizontal_cells * num_vertical_cells) {
            gridnumber = 0;
            //increase_overlays();
          }
        }
      }
    }
    pthread_mutex_unlock(&mutex);
    }
}

void pointsampler_mutate_with_pixel(path_t *curr, path_t *tent, float i, float j)
{
  path_init(tent, tent->index, tent->sensor.camid);
  path_set_pixel(tent, i, j);
  sampler_create_path(tent);
}

void pointsampler_enable_fake_random(pointsampler_t *s)
{
  const int tid = common_get_threadid();
  s->rand[tid].enabled = 1;
}

void pointsampler_disable_fake_random(pointsampler_t *s)
{
  const int tid = common_get_threadid();
  s->rand[tid].enabled = 0;
}

void pointsampler_set_fake_random(pointsampler_t *s, int dim, float rand)
{
  const int tid = common_get_threadid();
  s->rand[tid].rand[dim] = rand;
}

void pointsampler_prepare_frame(pointsampler_t *s)
{
  // would wrap around int limit and run out of bits for radical inverse?  stop
  // to re-init the random bit permutations (we only pass a seed, will be
  // randomised internally):
  if(rt.threads->end >> 32 > s->reinit)
    halton_init_random(&s->h, rt.anim_frame + ++s->reinit);
}

void pointsampler_stop_learning(pointsampler_t *s) { }
