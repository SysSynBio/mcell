#include "volume_output.h"
#include "mcell_structs.h"
#include "sched_util.h"
#include "vol_util.h"
#include "strfunc.h"

#include <errno.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

static void ioerror(FILE *err_file,
                    char const *fmt,
                    ...);

static int produce_item_header(FILE *err_file,
                               FILE *out_file,
                               struct volume_output_item *vo);

static int produce_mol_counts(struct volume *wrld,
                              FILE *out_file,
                              struct volume_output_item *vo);

static int find_species_in_array(struct species **mols,
                                 int num_mols,
                                 struct species *ptr);

static int reschedule_volume_output_item(struct volume *wrld,
                                         struct volume_output_item *vo);

/*
 * Output a block of volume data as requested by the 'vo' object.
 */
int update_volume_output(struct volume *wrld, struct volume_output_item *vo)
{
  int failure = 0;
  char *filename;
  no_printf("Updating volume output at time %lld of %lld\n", wrld->it_time, wrld->iterations);

  /* build the filename */
  filename = alloc_sprintf("%s.%lld.dat", vo->filename_prefix, wrld->it_time);
  if (filename == NULL) {
    fprintf(wrld->err_file,"File '%s', Line %ld: Out of memory while formatting volume output filename '%s'\n", __FILE__, (long)__LINE__, vo->filename_prefix);
    return 1;
  }

  /* Output the volume item */
  failure = output_volume_output_item(wrld, filename, vo);
  free(filename);

  /* Reschedule this volume item, if appropriate */
  if (! failure)
    failure = reschedule_volume_output_item(wrld, vo);

  /* Should we return failure if we can't create the file?  Doing so will bring
   * down the entire sim...
   */
  /* return failure; */
  return 0;
}

/*
 * Produce the output for a volume item.
 */
int output_volume_output_item(struct volume *wrld,
                              char const *filename,
                              struct volume_output_item *vo)
{
  FILE *f = fopen(filename, "w");
  if (f == NULL) {
    ioerror(wrld->err_file, "Couldn't open volume output file '%s'", filename);
    return 1;
  }

  if (produce_item_header(wrld->err_file, f, vo))
    goto failure;

  if (produce_mol_counts(wrld, f, vo))
    goto failure;

  fclose(f);
  return 0;

failure:
  fclose(f);
  return 1;
}

/*
 * Write the molecule counts to the file.
 *
 * XXX: Update this code to be smarter about the tradeoff between large slab
 * size (= large mem usage) and small slab size (= worse cache usage).
 */
static int produce_mol_counts(struct volume *wrld,
                              FILE *out_file,
                              struct volume_output_item *vo)
{
  struct volume_molecule *curmol;
  int *counters, *countersptr;
  int k, u, v;
  double z = vo->location.z, y = vo->location.y, x = vo->location.x;
  double x_lim = x + vo->voxel_size.x * (double) vo->nvoxels_x;
  double y_lim = y + vo->voxel_size.y * (double) vo->nvoxels_y;
  struct subvolume *cur_partition_z;
  double z_lim_part;

  /* Allocate memory for counters. */
  counters = (int *) malloc(sizeof(int) * vo->nvoxels_x * vo->nvoxels_y);

  cur_partition_z = find_subvolume(& vo->location, NULL);
  if (cur_partition_z == NULL)
  {
    fprintf(wrld->err_file, "While counting at [%g, %g, %g]: point isn't within a partition.", x, y, z);
    free(counters);
    return 1;
  }

  z_lim_part = wrld->z_fineparts[cur_partition_z->urb.z];

  /* For each slab: */
  double r_voxsz_x = 1.0 / vo->voxel_size.x;
  double r_voxsz_y = 1.0 / vo->voxel_size.y;
  for (k = 0; k < vo->nvoxels_z; ++k)
  {
    double z_lim_slab = z + vo->voxel_size.z;
    struct subvolume *cur_partition_y = cur_partition_z;

    /* reset counters for this slab */
    memset(counters, 0, sizeof(int) * vo->nvoxels_x * vo->nvoxels_y);

    /* Loop over relevant partitions */
keep_counting:
    cur_partition_y = cur_partition_z;
    while (cur_partition_y != NULL  &&  wrld->y_fineparts[cur_partition_y->llf.y] < y_lim)
    {
      struct subvolume *cur_partition = cur_partition_y;
      while (cur_partition != NULL  &&  wrld->x_fineparts[cur_partition->llf.x] < x_lim)
      {
        /* Count molecules */
        for (curmol = cur_partition->mol_head;
             curmol != NULL;
             curmol = curmol->next_v)
        {
          /* Skip defunct molecules */
          if (curmol->properties == NULL)
            continue;

          /* Skip molecules not in our slab */
          if (curmol->pos.z < z  ||  curmol->pos.z >= z_lim_slab)
            continue;

          /* Skip molecules outside our domain */
          if (curmol->pos.x < x  ||  curmol->pos.x >= x_lim  ||
              curmol->pos.x < y  ||  curmol->pos.y >= y_lim)
            continue;

          /* See if we're interested in this molecule */
          if (vo->num_molecules == 1) {
            if (*vo->molecules != curmol->properties) continue;
          }
          else {
            if ((find_species_in_array(vo->molecules, vo->num_molecules, curmol->properties)) == -1)
              continue;
          }

          /* We've got a winner!  Add one to the appropriate voxel. */
          ++ counters[((int) floor((curmol->pos.y - y) * r_voxsz_y)) * vo->nvoxels_y +
                       (int) floor((curmol->pos.x - x) * r_voxsz_x)];
        }

        /* Advance to next x-partition */
        cur_partition = traverse_subvol(cur_partition, NULL, X_POS);
      }

      /* Advance to next y-partition */
      cur_partition_y = traverse_subvol(cur_partition_y, NULL, Y_POS);
    }

    /* If the slab crosses a Z boundary, keep on truckin' */
    if (z_lim_slab >= z_lim_part)
    {
      /* If we can get to the next partition, don't update slab and don't
       * spill!
       */
      cur_partition_z = traverse_subvol(cur_partition_z, NULL, Z_POS);

      if (cur_partition_z != NULL) {
        z_lim_part = wrld->z_fineparts[cur_partition_z->urb.z];
        goto keep_counting;
      }
    }
    else
    {
      z = z_lim_slab;
      z_lim_slab += vo->voxel_size.z;
    }

    /* Spill our counts */
    countersptr = counters;
    for (u = 0; u < vo->nvoxels_y; ++u) {
      for (v = 0; v < vo->nvoxels_x; ++v)
        fprintf(out_file, "%d ", *countersptr++);
      fprintf(out_file, "\n");
    }

    /* Extra newline to put visual separation between slabs */
    fprintf(out_file, "\n");
  }

  free(counters);
  return 0;
}

/*
 * Binary search for a pointer in an array of pointers.
 */
static int find_species_in_array(struct species **mols,
                                 int num_mols,
                                 struct species *ptr)
{
  int lo = 0, hi = num_mols;
  while (hi - lo > 1)
  {
    int mid = (hi + lo) / 2;
    if (mols[mid] > ptr) hi = mid;
    else if (mols[mid] < ptr) lo = mid;
    else return mid;
  }

  if (mols[lo] == ptr) return lo;
  else return -1;
}

/*
 * Write the item header to the file.
 */
static int produce_item_header(FILE *err_file,
                               FILE *out_file,
                               struct volume_output_item *vo)
{
  if (fprintf(out_file, "# nx=%d ny=%d nz=%d time=%g\n",
              vo->nvoxels_x, vo->nvoxels_y, vo->nvoxels_z, vo->t) < 0)
  {
    ioerror(err_file, "Couldn't write header of volume output file");
    return 1;
  }

  return 0;
}

/*
 * Reschedule a volume output item, if necessary.
 */
static int reschedule_volume_output_item(struct volume *wrld,
                                         struct volume_output_item *vo)
{
  /* Find the next time */
  if (vo->timer_type == OUTPUT_BY_STEP)
    vo->t += (vo->step_time / wrld->time_unit);
  else
  {
    double time_scale = 0.0;

    /* Check if we're done */
    if (vo->next_time == vo->times + vo->num_times)
    {
      free(vo->filename_prefix);
      free(vo->molecules);
      free(vo->times);
      free(vo);
      return 0;
    }

    /* Compute the next time and advance the next_time ptr */
    if (vo->timer_type == OUTPUT_BY_ITERATION_LIST)
      time_scale = 1.0;
    else
      time_scale = 1.0 / wrld->time_unit;
    vo->t = (*vo->next_time++) * time_scale;
  }

  /* Add to the schedule */
  if (schedule_add(wrld->volume_output_scheduler, vo)) {
    fprintf(wrld->err_file,"File %s, Line %ld: Out of memory while setting up volume output.\n", __FILE__, (long)__LINE__);
    return 1;
  }

  return 0;
}

/*
 * Format an appropriate error message upon I/O failure.
 */
static void ioerror(FILE *err_file, char const *fmt, ...)
{
  va_list args;
  char buffer[1024];
  if (strerror_r(errno, buffer, sizeof(buffer)) == 0)
  {
    va_start(args, fmt);
    vfprintf(err_file, fmt, args);
    fprintf(err_file, ": %s\n", buffer);
    va_end(args);
  }
  else
  {
    va_start(args, fmt);
    vfprintf(err_file, fmt, args);
    fprintf(err_file, "\n");
    va_end(args);
  }
}

