
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <netdb.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <sys/resource.h>

#ifdef KELP
#include <kelp.h>
#endif

#ifndef RAN4_H
#include "ran4.h"
#endif

#include "rng.h"
#include "mcell_structs.h"
#include "strfunc.h"
#include "vector.h"
#include "sym_table.h"
#include "count_util.h"
#include "vol_util.h"
#include "wall_util.h"
#include "grid_util.h"
#include "viz_output.h"
#include "react.h"
#include "react_output.h"
#include "util.h"
#include "chkpt.h"
#include "mdlparse_util.h"
#include "init.h"

#ifdef DEBUG
#define no_printf printf
#endif

extern struct volume *world;

#define MICROSEC_PER_YEAR 365.25*86400.0*1e6

/**
 * Prints out author and grant credits.
 * When invoked with the -info option, this function prints to log_file
 * a header like:

 <pre>
  MCell (tm) Version 2.50x  06/05/2001  Running on dalton.salk.edu

  Copyright (C) 1997,1998,1999 by The Salk Institute & Cornell University
  Co-authored by Thomas M. Bartol Jr. & Joel R. Stiles
  Acknowledgements:
    The authors thank Edwin E. Salpeter for input on theory and
    algorithm design, Terrence J. Sejnowski for development input
    and support (NSF Grant IBN-9603611), and Miriam M. Salpeter
    for fostering quantitative experimental applications.
    Additional support from NIH Grant K08NS01776 (Joel R. Stiles).
 </pre>

 */
void init_credits(void)
{
  FILE *log_file;
  unsigned int seed;
  time_t the_time;
  char *institute[2],*author[2];
  int i;
  struct ran4_state rng;
 
  log_file=world->log_file;
  time(&the_time);
  seed=(unsigned int)the_time;
  ran4_init(&rng,seed);
  for (i=0;i<100;i++) ran4_uint32(&rng);
  if (ran4_dbl32(&rng)<0.5) {
    institute[0]=my_strdup("The Salk Institute");
    institute[1]=my_strdup("& Cornell University");
    author[0]=my_strdup("Thomas M. Bartol Jr.");
    author[1]=my_strdup("& Joel R. Stiles");
  }
  else {
    institute[0]=my_strdup("Cornell University");
    institute[1]=my_strdup("& The Salk Institute");
    author[0]=my_strdup("Joel R. Stiles");
    author[1]=my_strdup("& Thomas M. Bartol Jr.");
  }
  if((institute[0] == NULL) || (institute[1] == NULL) ||
     (author[0] == NULL) || (author[1] == NULL))
  {
     fprintf(stderr, "Out of memory.\n");
     return;
  }	
  
    fprintf(log_file,"  Copyright (C) 1997,1998,1999 by %s %s\n",institute[0],institute[1]);
    fprintf(log_file,"  Co-authored by %s %s\n",author[0],author[1]);
    fprintf(log_file,"  Acknowledgements:\n");
    fprintf(log_file,"    The authors thank Edwin E. Salpeter for input on theory and\n");
    fprintf(log_file,"    algorithm design, Terrence J. Sejnowski for development input\n");
    fprintf(log_file,"    and support (NSF Grant IBN-9603611), and Miriam M. Salpeter\n");
    fprintf(log_file,"    for fostering quantitative experimental applications.\n");
    fprintf(log_file,"    Additional support from NIH Grant K08NS01776 (Joel R. Stiles).\n\n");

  return;
}


/* Sets default notification values */
int init_notifications()
{
  world->notify = (struct notifications*)malloc(sizeof(struct notifications));
  if (world->notify==NULL) return 1;
 
  /* Notifications */
  world->notify->progress_report = NOTIFY_FULL;
  world->notify->diffusion_constants = NOTIFY_BRIEF;
  world->notify->reaction_probabilities = NOTIFY_FULL;
  world->notify->reaction_prob_notify = 0.0;
  world->notify->partition_location = NOTIFY_NONE;
  world->notify->box_triangulation = NOTIFY_NONE;
  world->notify->custom_iterations = NOTIFY_FULL;
  world->notify->custom_iteration_value = 0;  /* Ignored unless NOTIFY_CUSTOM set */
  world->notify->release_events = NOTIFY_FULL;
  world->notify->file_writes = NOTIFY_NONE;
  world->notify->final_summary = NOTIFY_FULL;
  /* Warnings */
  world->notify->neg_diffusion = WARN_WARN;
  world->notify->neg_reaction = WARN_WARN;
  world->notify->high_reaction_prob = WARN_COPE;
  world->notify->reaction_prob_warn = 1.0;
  world->notify->close_partitions = WARN_WARN;
  world->notify->degenerate_polys = WARN_WARN;
  world->notify->overwritten_file = WARN_COPE;
  world->notify->short_lifetime = WARN_WARN;
  world->notify->short_lifetime_value = 50;
  world->notify->missed_reactions = WARN_WARN;
  world->notify->missed_reaction_value = 0.001;
  world->notify->missed_surf_orient = WARN_ERROR;
  world->notify->useless_vol_orient = WARN_WARN;
  
  if (world->log_freq!=-1) /* User set this */
  {
    world->notify->custom_iterations = NOTIFY_CUSTOM;
    world->notify->custom_iteration_value = world->log_freq;
  }

  return 0;
}



/**
 * Initializes the parameters required for the simulation (duh!).
 * This driver function sets up everything needed to get the simulation
 * up and running. It does the following (in that approximate order):
 * 	- Initializes variables to nice values
 * 	- Initialize all the compute nodes (init_nodes())
 * 	- Creates all the data structures
 * 	- Parse the MDL input file (mdlparse_init())
 * 	- Setup the geometry (init_geom())
 * \todo Need more info here.
 */
int init_sim(void)
{
  FILE *log_file;
  struct sym_table *gp;
  struct output_block *obp,*obpn;
  int i;
  int *intp;
  int reactants_3D_present = 0; /* flag to check whether there are 3D reactants
                             (participants in the reactions
                              between 3D molecules) in the simulation */

#ifdef USE_RAN4
#include "seed_array.h"
#endif 

  log_file=world->log_file;
#ifdef KELP
  if (world->procnum == 0) {
#endif
  fprintf(log_file,"MCell initializing simulation...\n");
  fflush(log_file);
#ifdef KELP
  }
#endif


  /* Initialize variables to reasonably safe values */

  world->curr_file=world->mdl_infile_name;

  /* by Erhan Gokcay 5/3/2002 ========================== */
  /* We can not initialize chkpt_infile anymore. It is set in mcell.c */
  /*  chkpt_infile=NULL; */
  /* =================================================== */

  world->chkpt_outfile=NULL;
  world->chkpt_iterations=0;
  world->chkpt_seq_num=0;

  /*world->chkpt_init=1; */  /* set in the main() */
  world->chkpt_flag=0;
  world->molecule_prefix_name=NULL;
  world->file_prefix_name=NULL;
  world->random_number_use=0;
  world->ray_voxel_tests=0;
  world->ray_polygon_tests=0;
  world->ray_polygon_colls=0;
  world->mol_mol_colls=0;
  world->diffusion_steps=0;
  world->sim_elapsed_time=0;
  world->chkpt_elapsed_real_time=0;
  world->chkpt_elapsed_real_time_start=0;
  world->chkpt_byte_order_mismatch = 0;
  world->it_time=0;
  world->elapsed_time=0;
  world->time_unit=0;
  world->time_step_max=0;
  world->start_time=0;
  world->current_real_time=0;
  world->current_start_real_time=0;
  world->effector_grid_density=10000;
  world->length_unit=1.0/sqrt(world->effector_grid_density);
  world->rx_radius_3d = 0;
  world->max_diffusion_step=0;
  world->radial_directions=16384;
  world->radial_subdivisions=1024;
  world->fully_random=0;
  world->num_directions=world->radial_directions;
  world->r_num_directions=1.0/world->num_directions;
  world->r_step=NULL;
  world->r_step_surface=NULL;
  world->d_step=NULL;
  world->place_waypoints_flag=0;
  world->releases_on_regions_flag=0;
  world->count_scheduler = NULL;
  world->storage_head = NULL;
  world->storage_allocator = NULL;
  world->x_partitions = NULL;
  world->y_partitions = NULL;
  world->z_partitions = NULL;
  world->x_fineparts = NULL;
  world->y_fineparts = NULL;
  world->z_fineparts = NULL;
  world->n_fineparts = 0;
  
  world->viz_output_flag = 0; 
  world->use_expanded_list=1;
  world->randomize_gmol_pos=0;
  world->vacancy_search_dist2=0;
  world->surface_reversibility=0;
  world->volume_reversibility=0;
  
  world->mcell_version = MCELL_VERSION;
  
  world->clamp_list = NULL;

  world->rng = malloc(sizeof(struct rng_state));
  if (world->rng==NULL)
  {
    fprintf(world->err_file,"Out of memory: failed to allocate random number generator\n");
    exit(EXIT_FAILURE);
  }
#ifdef USE_RAN4
  if (world->seed_seq < 1 || world->seed_seq > 3000) {
    fprintf(log_file,"MCell: error, random sequence number not in range 1 to 3000\n  Recompile without USE_RAN4 flag in rng.h to increase the range\n");
    return(1);
  }
  world->init_seed = seed_array[world->seed_seq-1];
  rng_init(world->rng,world->init_seed);
  fprintf(log_file,"MCell[%d]: random sequence: %d  seed: %d\n", world->procnum,world->seed_seq,world->init_seed);
  fflush(log_file);
#else
  if (world->seed_seq < 1 || world->seed_seq > INT_MAX) {
    fprintf(log_file,"MCell: error, random sequence number not in range 1 to 2^31-1\n");
    return(1);
  }
  rng_init(world->rng,world->seed_seq);
  fprintf(log_file,"MCell[%d]: random sequence %d\n",world->procnum,world->seed_seq);
  fflush(log_file);
#endif

  world->count_hashmask = COUNT_HASHMASK;
  world->count_hash = (struct counter**)malloc(sizeof(struct counter*)*(world->count_hashmask+1));
  if (world->count_hash == NULL)
  {
    fprintf(log_file,"MCell: Out of memory while creating counter hash table\n");
    exit(EXIT_FAILURE);
  }
  for (i=0;i<=world->count_hashmask;i++) world->count_hash[i] = NULL;
  
  world->pathway_requester = create_mem(sizeof(struct pathway_count_request),64);
  if (world->pathway_requester==NULL)
  {
    fprintf(world->err_file,"Out of memory: could not create space to pair reactions with count requests\n");
    exit(EXIT_FAILURE);
  }

  if((world->main_sym_table=init_symtab(SYM_HASHSIZE)) == NULL){
    fprintf(log_file,"MCell: initialization of symbol table failed\n");
    return(1);
  }
	

  if ((gp=store_sym("WORLD_OBJ",OBJ,world->main_sym_table))==NULL) {
    fprintf(log_file,"MCell: Out of memory while creating world root object\n");
    return(1);
  }
  world->root_object=(struct object *)gp->value;
  world->root_object->object_type=META_OBJ;
  world->root_object->last_name="";

  if ((gp=store_sym("WORLD_INSTANCE",OBJ,world->main_sym_table))==NULL) {
    fprintf(log_file,"MCell: Out of memory while creating world root instance\n");
    return(1);
  }
  world->root_instance=(struct object *)gp->value;
  world->root_instance->object_type=META_OBJ;
  world->root_instance->last_name="";

  if ((gp=store_sym("DEFAULT_RELEASE_PATTERN",RPAT,world->main_sym_table))
      ==NULL) {
    fprintf(log_file,"MCell: Out of memory while creating default release pattern");
    return(1);
  }
  world->default_release_pattern=(struct release_pattern *)gp->value;
  world->default_release_pattern->delay=0;
  world->default_release_pattern->release_interval=FOREVER;
  world->default_release_pattern->train_interval=FOREVER;
  world->default_release_pattern->train_duration=FOREVER;
  world->default_release_pattern->number_of_trains=1;
   
  if ((gp=store_sym("GENERIC_MOLECULE",MOL,world->main_sym_table))
      ==NULL) {
    fprintf(log_file,"MCell: Out of memory while creating generic molecule");
    return(1);
  }
  world->g_mol=(struct species *)gp->value;

  if ((gp=store_sym("GENERIC_SURFACE",MOL,world->main_sym_table))
      ==NULL) {
    fprintf(log_file,"MCell: Out of memory while creating generic surface");
    return(1);
  }
  world->g_surf=(struct species *)gp->value;
  world->g_surf->flags=IS_SURFACE;

  world->output_block_head=NULL;
  world->release_event_queue_head=NULL;
  world->tot_mols=0;
  world->viz_obj_head=NULL;
  world->viz_mode=DREAMM_V3_MODE;
  world->rk_mode_var=NULL;
  world->frame_data_head=NULL;

  if ((world->count_zero=(struct output_evaluator *)malloc
       (sizeof(struct output_evaluator)))==NULL) {
    fprintf(log_file,"MCell: Out of memory while creating zero counter\n");
    exit(EXIT_FAILURE);
  }
  if (!(intp=(int *)malloc(sizeof(int)))) {
    fprintf(log_file,"MCell: Out of memory while creating zero counter\n");
    exit(EXIT_FAILURE);
  }
  *intp=0;
  world->count_zero->next=NULL;
  world->count_zero->update_flag=0;
  world->count_zero->reset_flag=0;
  world->count_zero->index_type=TIME_STAMP_VAL;
  world->count_zero->n_data=1;
  world->count_zero->data_type=INT;
  world->count_zero->temp_data=(void *)intp;
  world->count_zero->final_data=(void *)intp;
  world->count_zero->operand1=NULL;
  world->count_zero->operand2=NULL;
  world->count_zero->oper='\0';
  
  world->releaser = create_scheduler(1.0,100.0,100,0.0);
  if(world->releaser == NULL){
	fprintf(stderr, "Out of memory while creating releaser.\n");
        exit(EXIT_FAILURE);
  }

  /* Parse the MDL file: */
  no_printf("Node %d parsing MDL file %s\n",world->procnum,world->mdl_infile_name);
  fflush(stderr);
  if (mdlparse_init(world)) {
    fprintf(log_file,"MCell: error parsing file: %s\n",world->curr_file);
    return(1);
  }
  no_printf("Done parsing MDL file: %s\n",world->mdl_infile_name);
  fflush(stderr);

  /* Set up the array of species */
  if (init_species())
  {
    fprintf(log_file,"MCell: error initializing species\n");
    return(1);
  }
  no_printf("Done setting up species.\n");
  

  /* Visualize all molecules if asked in "mdl" file */
  if((world->viz_mode == DREAMM_V3_MODE) || (world->viz_mode == DREAMM_V3_GROUPED_MODE))
  {
    if((world->viz_output_flag & VIZ_ALL_MOLECULES) != 0) {
       struct species *sp;
  
      for(i = 0; i < world->n_species; i++)
      {
         sp = world->species_list[i];
         if((sp->flags & IS_SURFACE) != 0) continue;
         if(strcmp(sp->sym->name, "GENERIC_MOLECULE") == 0) continue;  

         /* set viz_state to INCLUDE_OBJ for the molecule we want to visualize
             but will not assign state value */
         sp->viz_state = INCLUDE_OBJ;

      } 

    }
  }


 /* If there are no 3D molecules-reactants in the simulation
    set up the"use_expanded_list" flag to zero. */
  for(i = 0; i < world->n_species; i++)
  {
        struct	species *sp = world->species_list[i];
        if((sp->flags & CAN_MOLMOL) != 0){
		reactants_3D_present = 1;
                break;
        }
  }
  if(reactants_3D_present == 0){
	world->use_expanded_list = 0;
  }

/* Instantiation Pass #1: Initialize the geometry */
  if (init_geom()) {
    fprintf(log_file,"MCell: error initializing geometry\n");
    return(1);
  }
  
  no_printf("Done setting up geometry.\n");
  
/* Instantiation Pass #2: Partition geometry */
  if (init_partitions()) {
    fprintf(log_file,"MCell: error initializing partitions.\n");
    return(1);
  }
  
  if (distribute_world()) {
    fprintf(log_file,"MCell: error moving geometry to partitions\n");
    return(1);
  }
  
  if (sharpen_world()) {
    fprintf(log_file,"MCell: error adding edges to geometry\n");
    return(1);
  }

/* Instantiation Pass #3: Initialize regions */
  if (init_regions()) {
    fprintf(log_file,"MCell: error initializing object regions\n");
    return(1);
  }
  
  if (check_region_counters()) {
    fprintf(log_file,"MCell: error in region count statement\n");
    return(1);
  }

  if (world->place_waypoints_flag) {
    if (place_waypoints()) {
      fprintf(log_file,"MCell: error storing waypoints\n");
      return(1);
    }
  }
  
  if (init_effectors())
  {
    fprintf(world->err_file,"Error initializing effectors on regions.\n");
    return 1;
  }
  
  if (world->releases_on_regions_flag)
  {
    if (init_releases())
    {
      fprintf(log_file,"MCell: error intializing releases on regions\n");
      return 1;
    }
  }

  /* Decompose the space into subvolumes */
/*
  if (decompose_volume(volume,wall_head)) {
    fprintf(log_file,"MCell: error decomposing volume\n");
    return(1);
  }
  resolve_contiguity(cmprt_head);
  if (init_effector_table(wall_head)) {
    fprintf(log_file,"MCell: error initializing effector table\n");
    return(1);
  }
  fflush(log_file);
  if (init_rx_table(rx_head)) {
    fprintf(log_file,"MCell: error initializing rx table\n");
    return(1);
  }
  if (init_release_event_table(release_event_queue_head)) {
    fprintf(log_file,"MCell: error initializing release event table\n");
    return(1);
  }
*/  

  if (world->chkpt_infile) {
    if ((world->chkpt_infs=fopen(world->chkpt_infile,"rb"))==NULL) {
      world->chkpt_seq_num=1;
    }
    else {
      fprintf(log_file,"MCell: reading from checkpoint file %s\n",world->chkpt_infile);
      if(read_chkpt(world->chkpt_infs)) {
	fprintf(log_file,"MCell: error reading from checkpoint file %s\n",world->chkpt_infile);
	return(1);
      }
      fclose(world->chkpt_infs);
    }
  }
  else {
    world->chkpt_seq_num=1;
  }

   /**
   *Initialize the frame data list for the visualization 
   *and reaction output.
   **/
  init_frame_data_list(world->frame_data_head); 
/*
  init_reaction_list(reaction_data_head);
*/


  world->count_scheduler = create_scheduler(1.0,100.0,100,world->current_start_real_time/world->time_unit);
  if(world->count_scheduler == NULL){
	fprintf(stderr,"Out of memory while creating count_scheduler.\n");
        exit(EXIT_FAILURE);
  }

/* Schedule the reaction data output events */
  obp = world->output_block_head;
  while(obp != NULL)
  {
    obpn = obp->next;
    if(obp->timer_type == OUTPUT_BY_STEP)
    {
       if(world->chkpt_seq_num > 1)
       {
          while((obp->t < world->iterations + 1) && (obp->t <= world->count_scheduler->now)){
              obp->t += obp->step_time/world->time_unit;

          }
       }      
    }else if(obp->timer_type != OUTPUT_BY_STEP && obp->curr_time_ptr == NULL){
        obp->curr_time_ptr = obp->time_list_head;
        if(world->chkpt_seq_num == 1){
          if(obp->timer_type == OUTPUT_BY_ITERATION_LIST){
             obp->t = obp->curr_time_ptr->value;
          }else{
             obp->t = obp->curr_time_ptr->value/world->time_unit;
          }
        }else{
          if(obp->timer_type == OUTPUT_BY_ITERATION_LIST){
             while((obp->t < world->iterations + 1) && (obp->t <= world->count_scheduler->now)){
                 obp->curr_time_ptr = obp->curr_time_ptr->next;
                 if(obp->curr_time_ptr == NULL) break;
                 obp->t = obp->curr_time_ptr->value;
              }
          }else{
             while((obp->t < world->iterations + 1) && (obp->t <= world->count_scheduler->now)){
                 obp->curr_time_ptr = obp->curr_time_ptr->next;
                 if(obp->curr_time_ptr == NULL) break;
                 obp->t = obp->curr_time_ptr->value/world->time_unit;
              }
          }
        }
    }
    if (schedule_add(world->count_scheduler , obp))
    {
      	fprintf(stderr,"Out of memory: trying to save intermediate results\n");
      	int i = emergency_output();
      	fprintf(stderr,"Fatal error: out of memory while scheduling output for count statements.\nAttempt to write intermediate results had %d errors\n", i);
      	exit(EXIT_FAILURE);
    }
    obp = obpn;
  }

  no_printf("Done initializing simulation\n");
  fflush(log_file);
  return(0);
}


/********************************************************************
init_species: 
   Initializes array of molecules types to the default properties values.

*********************************************************************/
int init_species(void)
{
  int i;
  int count = 0;
  struct sym_table *gp;
  struct species *s;
  double speed;
  
  world->speed_limit = 0;
  
  for (i=0;i<SYM_HASHSIZE;i++)
  {
    for (gp = world->main_sym_table[i] ; gp != NULL ; gp = gp->next)
    {    
      if (gp->sym_type==MOL) count++;
    }
  }
  
  world->n_species = count;
  if((world->species_list = (struct species**)malloc(sizeof(struct species*)*world->n_species)) == NULL)
  {
	fprintf(stderr, "Out of memory during species initialization.\n");
        exit(EXIT_FAILURE);
  }
  count = 0;
  for (i=0;i<SYM_HASHSIZE;i++)
  {
    for (gp = world->main_sym_table[i] ; gp != NULL ; gp = gp->next)
    {    
      if (gp->sym_type==MOL)
      {
        s = (struct species*) gp->value;
        world->species_list[count] = s;
        world->species_list[count]->species_id = count;
        world->species_list[count]->chkpt_species_id = UINT_MAX;
/*        world->species_list[count]->hashval &= world->rx_hashsize-1; */
        world->species_list[count]->radius = EPS_C;
        world->species_list[count]->population = 0;
	world->species_list[count]->n_deceased = 0;
	world->species_list[count]->cum_lifetime = 0;
        if(world->species_list[count]->viz_state < 0){
        	world->species_list[count]->viz_state = EXCLUDE_OBJ;
        }
        if ( (s->flags & NOT_FREE) == 0 )
        {
          speed = 6.0*s->space_step/sqrt(MY_PI);
          if (speed > world->speed_limit) world->speed_limit = speed;
        }
        count++;
      }
    }
  }
 
   
  return 0;
}


/* This is just a placeholder for now. */
int init_partitions(void)
{
  int i,j,k,h;
  struct vector3 part_lo,part_hi;
  struct subvolume *sv;
  struct storage *shared_mem;
  
  if (world->bb_min.x >= world->bb_max.x)
  {
    part_lo.x = -10.0;
    part_hi.x = 10.0;
  }
  else
  {
    part_lo.x = world->bb_min.x;
    part_hi.x = world->bb_max.x;
  }
  if (world->bb_min.y >= world->bb_max.y)
  {
    part_lo.y = -10.0;
    part_hi.y = 10.0;
  }
  else
  {
    part_lo.y = world->bb_min.y;
    part_hi.y = world->bb_max.y;
  }
  if (world->bb_min.z >= world->bb_max.z)
  {
    part_lo.z = -10.0;
    part_hi.z = 10.0;
  }
  else
  {
    part_lo.z = world->bb_min.z;
    part_hi.z = world->bb_max.z;
  }

/* Cheating for now. */  
  part_lo.x = -1.1;
  part_hi.x = 1.01;
  part_lo.y = -1.1;
  part_hi.y = 1.01;
  part_lo.z = -1.1;
  part_hi.z = 1.01;

  part_lo.x *= (1.0+EPS_C)/world->length_unit;
  part_lo.y *= (1.0+EPS_C)/world->length_unit;
  part_lo.z *= (1.0+EPS_C)/world->length_unit;
  part_hi.x *= (1.0+EPS_C)/world->length_unit;
  part_hi.y *= (1.0+EPS_C)/world->length_unit;
  part_hi.z *= (1.0+EPS_C)/world->length_unit;
  
  
#if 0
  world->nx_parts = 11;

  if (world->nx_parts < 4) world->nx_parts = 4;
  
  world->ny_parts = world->nz_parts = world->nx_parts;
  
  world->x_partitions = (double*)malloc(sizeof(double)*world->nx_parts);
  world->y_partitions = (double*)malloc(sizeof(double)*world->ny_parts);
  world->z_partitions = (double*)malloc(sizeof(double)*world->nz_parts);

  world->x_partitions[0] = - GIGANTIC;
  world->x_partitions[1] = part_lo.x;
  world->x_partitions[world->nx_parts-2] = part_hi.x;
  world->x_partitions[world->nx_parts-1] = GIGANTIC;

  world->y_partitions[0] = - GIGANTIC;
  world->y_partitions[1] = part_lo.y;
  world->y_partitions[world->ny_parts-2] = part_hi.y;
  world->y_partitions[world->ny_parts-1] = GIGANTIC;

  world->z_partitions[0] = - GIGANTIC;
  world->z_partitions[1] = part_lo.z;
  world->z_partitions[world->nz_parts-2] = part_hi.z;
  world->z_partitions[world->nz_parts-1] = GIGANTIC;

  for (i=2;i<world->n_parts-2;i++)
  {
    frac = ((double)(i-1)) / ((double)(world->n_parts-3));
    world->x_partitions[i] = part_lo.x + frac*(part_hi.x - part_lo.x);
    world->y_partitions[i] = part_lo.y + frac*(part_hi.y - part_lo.y);
    world->z_partitions[i] = part_lo.z + frac*(part_hi.z - part_lo.z);
  }
    
  world->n_fineparts = world->n_parts;
  world->x_fineparts = world->x_partitions;
  world->y_fineparts = world->y_partitions;
  world->z_fineparts = world->z_partitions;
#endif

  if (set_partitions()) return 1;
  
  world->n_waypoints = 1;
  if((world->waypoints = (struct waypoint*)malloc(sizeof(struct waypoint*)*world->n_waypoints)) == NULL)
  {
    fprintf(stderr,"Fatal error: out of memory while initializing partitions.\n");
    exit(EXIT_FAILURE);
  }
  
  if((shared_mem = (struct storage*)malloc(sizeof(struct storage))) == NULL)
  {
    fprintf(stderr,"Fatal error: out of memory while initializing partitions.\n");
    exit(EXIT_FAILURE);
  }
  
  shared_mem->list = create_mem(sizeof(struct wall_list),128);
  shared_mem->mol  = create_mem(sizeof(struct molecule),128);
  shared_mem->gmol  = create_mem(sizeof(struct grid_molecule),128);
  shared_mem->face = create_mem(sizeof(struct wall),128);
  shared_mem->join = create_mem(sizeof(struct edge),128);
  shared_mem->tree = create_mem(sizeof(struct vertex_tree),128);
  shared_mem->effs = create_mem(sizeof(struct surface_grid),128);
  shared_mem->coll = create_mem(sizeof(struct collision),128);
  shared_mem->regl = create_mem(sizeof(struct region_list),128);
  shared_mem->exdv = create_mem(sizeof(struct exd_vertex),64);
  
  if (shared_mem->list==NULL ||
      shared_mem->mol==NULL  || shared_mem->gmol==NULL ||
      shared_mem->face==NULL || shared_mem->join==NULL ||
      shared_mem->tree==NULL || shared_mem->effs==NULL ||
      shared_mem->coll==NULL || shared_mem->regl==NULL ||
      shared_mem->exdv==NULL)
  {
    fprintf(stderr,"Fatal error: out of memory while initializing partitions.\n");
    exit(EXIT_FAILURE);
  }
  
  shared_mem->wall_head = NULL;
  shared_mem->wall_count = 0;
  shared_mem->vert_head = NULL;
  shared_mem->vert_count = 0;
 
 if(world->chkpt_init)
 {
     if((shared_mem->timer = create_scheduler(1.0,100.0,100,0.0)) == NULL){
	fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        int i = emergency_output();
	fprintf(stderr,"Fatal error: out of memory while partition initialization.\nAttempt to write intermediate results had %d errors.\n", i);
        exit(EXIT_FAILURE);
     }
     shared_mem->current_time = 0.0;
  }
  
  if (world->time_step_max==0.0) shared_mem->max_timestep = MICROSEC_PER_YEAR;
  else
  {
    if (world->time_step_max < world->time_unit) shared_mem->max_timestep = 1.0;
    else shared_mem->max_timestep = world->time_step_max/world->time_unit;
  }
  
  if((world->storage_allocator = create_mem(sizeof(struct storage_list),10)) == NULL)
  {
	fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        int i = emergency_output();
	fprintf(stderr,"Fatal error: out of memory while partition initialization.\nAttempt to write intermediate results had %d errors.\n", i);
        exit(EXIT_FAILURE);
  }
  if((world->storage_head = (struct storage_list*)mem_get(world->storage_allocator)) == NULL) {
	fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        int i = emergency_output();
	fprintf(stderr,"Fatal error: out of memory while partition initialization.\nAttempt to write intermediate results had %d errors.\n", i);
        exit(EXIT_FAILURE);
  }
  world->storage_head->next = NULL;
  world->storage_head->store = shared_mem;
  
  world->n_subvols = (world->nz_parts-1) * (world->ny_parts-1) * (world->nx_parts-1);
  printf("Creating %d subvolumes (%d,%d,%d per axis)\n",world->n_subvols,world->nx_parts-1,world->ny_parts-1,world->nz_parts-1);
  if((world->subvol = (struct subvolume*)malloc(sizeof(struct subvolume)*world->n_subvols)) == NULL){
	fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        int i = emergency_output();
	fprintf(stderr,"Fatal error: out of memory while partition initialization.\nAttempt to write intermediate results had %d errors.\n", i);
        exit(EXIT_FAILURE);
  }
  for (i=0;i<world->nx_parts-1;i++)
  for (j=0;j<world->ny_parts-1;j++)
  for (k=0;k<world->nz_parts-1;k++)
  {
    h = k + (world->nz_parts-1)*(j + (world->ny_parts-1)*i);
    sv = & (world->subvol[ h ]);
    sv->wall_head = NULL;
    sv->wall_tail = NULL;
    sv->wall_count = 0;
    sv->mol_head = NULL;
    sv->mol_count = 0;
    
    sv->index = h;
    
    sv->llf.x = bisect_near( world->x_fineparts , world->n_fineparts , world->x_partitions[i] );
    sv->llf.y = bisect_near( world->y_fineparts , world->n_fineparts , world->y_partitions[j] );
    sv->llf.z = bisect_near( world->z_fineparts , world->n_fineparts , world->z_partitions[k] );
    sv->urb.x = bisect_near( world->x_fineparts , world->n_fineparts , world->x_partitions[i+1] );
    sv->urb.y = bisect_near( world->y_fineparts , world->n_fineparts , world->y_partitions[j+1] );
    sv->urb.z = bisect_near( world->z_fineparts , world->n_fineparts , world->z_partitions[k+1] );
    
    sv->is_bsp = 0;

    if (i==0) sv->neighbor[X_NEG] = NULL;
    else sv->neighbor[X_NEG] = &(world->subvol[ h - (world->nz_parts-1)*(world->ny_parts-1) ]);
    
    if (i==world->nx_parts-2) sv->neighbor[X_POS] = NULL;
    else sv->neighbor[X_POS] = &(world->subvol[ h + (world->nz_parts-1)*(world->ny_parts-1) ]);
    
    if (j==0) sv->neighbor[Y_NEG] = NULL;
    else sv->neighbor[Y_NEG] = &(world->subvol[ h - (world->nz_parts-1) ]);
    
    if (j==world->ny_parts-2) sv->neighbor[Y_POS] = NULL;
    else sv->neighbor[Y_POS] = &(world->subvol[ h + (world->nz_parts-1) ]);
    
    if (k==0) sv->neighbor[Z_NEG] = NULL;
    else sv->neighbor[Z_NEG] = &(world->subvol[ h - 1 ]);
    
    if (k==world->nz_parts-2) sv->neighbor[Z_POS] = NULL;
    else sv->neighbor[Z_POS] = &(world->subvol[ h + 1 ]);
    
    sv->local_storage = shared_mem;
  }
  
  world->binning = 0;
  world->lookup = NULL;
  
  return 0;
}



/**
 * Initializes the geometry of the world.
 * Calls instance_obj() to instantiate all physical objects.
 * (Meta objects, box objects, polygon objects and release sites)
 * Populates viz_obj list vizp and lig_count_ref list lcrp.
 */
int init_geom(void)
{
  FILE *log_file;
  double tm[4][4];
  double vol_infinity;
  struct release_event_queue *req,*rqn;
  
  no_printf("Initializing physical objects\n");
  log_file=world->log_file;
  vol_infinity=sqrt(DBL_MAX)/4;
  world->bb_min.x=vol_infinity;
  world->bb_min.y=vol_infinity;
  world->bb_min.z=vol_infinity;
  world->bb_max.x=-vol_infinity;
  world->bb_max.y=-vol_infinity;
  world->bb_max.z=-vol_infinity;
  init_matrix(tm);
  
  compute_bb(world->root_instance,tm,NULL);
  if (world->bb_min.x==vol_infinity 
      && world->bb_min.y==vol_infinity
      && world->bb_min.z==vol_infinity
      && world->bb_max.x==-vol_infinity
      && world->bb_max.y==-vol_infinity
      && world->bb_max.z==-vol_infinity) {
    world->bb_min.x=0;
    world->bb_min.y=0;
    world->bb_min.z=0;
    world->bb_max.x=0;
    world->bb_max.y=0;
    world->bb_max.z=0;
  }
  if (world->procnum == 0) {
    fprintf(log_file,"MCell: world bounding box in microns =\n");
    fprintf(log_file,"         [ %.9g %.9g %.9g ] [ %.9g %.9g %.9g ]\n",
      world->bb_min.x*world->length_unit,world->bb_min.y*world->length_unit,
      world->bb_min.z*world->length_unit,world->bb_max.x*world->length_unit,
      world->bb_max.y*world->length_unit,world->bb_max.z*world->length_unit);
  }

  world->n_walls=world->root_instance->n_walls;
  world->n_verts=world->root_instance->n_verts;
  no_printf("World object contains %d walls and %d vertices\n",
    world->n_walls,world->n_verts);
  
  if (instance_obj(world->root_instance,tm,NULL,NULL,NULL)) {
    return(1);
  }
  
/* Stick the queue of release events in a scheduler */
  req = world->release_event_queue_head;  
  while(req != NULL)
  {
    rqn = req->next;
    if (schedule_add(world->releaser , req)){ 
	fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        int i = emergency_output();
	fprintf(stderr,"Fatal error: out of memory while geometry initialization.\nAttempt to write intermediate results had %d errors.\n", i);
        exit(EXIT_FAILURE);
    }
    req = rqn;
  }

  return(0);
}



/**
 * Instantiates all physical objects.
 * This function is recursively called on the tree object objp until
 * all the objects in the data structure have been instantiated.
 * <br>
 * This function actually calls instance_release_site() and
 * instance_polygon_object() to handle the actual instantiation of
 * those objects.
 */
int instance_obj(struct object *objp, double (*im)[4], struct viz_obj *vizp, struct lig_count_ref *lcrp, char *sub_name)
{
  FILE *log_file;
  struct object *child_objp;
  double tm[4][4];
  unsigned short l,m,n;
  char *tmp_name;
  int i;

  log_file=world->log_file;
  l=4;
  m=4;
  n=4;
  mult_matrix(objp->t_matrix,im,tm,l,m,n);
  if (vizp==NULL) {
    vizp=objp->viz_obj;
  }
  if (lcrp==NULL) {
    lcrp=objp->lig_count_ref;
  }

  if (sub_name!=NULL) { 
    if (strcmp(sub_name,"")==0) {
      tmp_name=my_strdup("");
      if(tmp_name == NULL){
	fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        i = emergency_output();
	fprintf(stderr,"Fatal error: out of memory while object instantiation.\nAttempt to write intermediate results had %d errors.\n", i);
        exit(EXIT_FAILURE);
      }else{}
    }
    else {
      tmp_name=my_strcat(sub_name,"."); 
      if(tmp_name == NULL){
	fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        i = emergency_output();
	fprintf(stderr,"Fatal error: out of memory while object instantiation.\nAttempt to write intermediate results had %d errors.\n", i);
        exit(EXIT_FAILURE);
      }else{}
                  
    }
    sub_name=my_strcat(tmp_name,objp->last_name);    
    if(sub_name == NULL){
	fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        i = emergency_output();
	fprintf(stderr,"Fatal error: out of memory while object instantiation.\nAttempt to write intermediate results had %d errors.\n", i);
        exit(EXIT_FAILURE);
    }else{}
    free((void *)tmp_name);
  }
  else {
    sub_name=my_strdup(objp->last_name);    
    if(sub_name == NULL){
	fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        i = emergency_output();
	fprintf(stderr,"Fatal error: out of memory while object instantiation.\nAttempt to write intermediate results had %d errors.\n", i);
        exit(EXIT_FAILURE);
    }else{}
  }
  
  switch (objp->object_type) {
  case META_OBJ:
    no_printf("Meta object %s instanced\n",sub_name);
    fflush(log_file);
    child_objp=objp->first_child;
    while (child_objp!=NULL) {
      if (instance_obj(child_objp,tm,vizp,lcrp,sub_name)) {
        return(1);
      }
      child_objp=child_objp->next;
    }
    break;
  case REL_SITE_OBJ:
    no_printf("Release site object %s instanced\n",sub_name);
    fflush(log_file);
    if (instance_release_site(objp,tm)) {
      return(1);
    }
    break;
  case BOX_OBJ:
    no_printf("Box object %s instanced\n",sub_name);
    fflush(log_file);
    if (instance_polygon_object(objp,tm,vizp,lcrp,sub_name)) {
      return(1);
    }
    break;
  case POLY_OBJ:
    no_printf("Polygon list object %s instanced\n",sub_name);
    fflush(log_file);
    if (instance_polygon_object(objp,tm,vizp,lcrp,sub_name)) {
      return(1);
    }
    break;
  }

  free((void *)sub_name);
  return(0);
}



/**
 * Instantiates a release site.
 * Creates a new release site from a template release site
 * as defined in the MDL file after applying the necessary
 * geometric transformations (rotation and translation).
 * Adds the rel
 */
int instance_release_site(struct object *objp, double (*im)[4])
{
  FILE *log_file;
  struct release_site_obj *rsop;
  struct release_event_queue *reqp;
  int i,j;

  log_file=world->log_file;
  rsop=(struct release_site_obj *)objp->contents;
  
  no_printf("Instancing release site object %s\n",objp->sym->name);
  fflush(log_file);
  reqp = (struct release_event_queue*)malloc(sizeof(struct release_event_queue));
  if (reqp==NULL)
  {
    fprintf(world->err_file,"Fatal error: out of memory while instantiating release site.\n");
    return 1;
  }

  reqp->release_site=rsop;
  reqp->event_time=rsop->pattern->delay;
  reqp->train_counter=0;
  reqp->train_high_time=rsop->pattern->delay;
  for (i=0;i<4;i++) for (j=0;j<4;j++) reqp->t_matrix[i][j]=im[i][j];
  reqp->next=world->release_event_queue_head;
  world->release_event_queue_head=reqp;
  

  if(rsop->pattern->train_duration > rsop->pattern->train_interval)
  {
    fprintf(world->err_file,"Error: Release pattern train duration is greater than train interval\n");
    return 1;
  } 
  
  no_printf("Done instancing release site object %s\n",objp->sym->name);

  fflush(log_file);
  return(0);
}



/**
 * Computes the bounding box for the entire simulation world.
 * Does things recursively in a manner similar to instance_obj().
 */
int compute_bb(struct object *objp, double (*im)[4], char *sub_name)
{
  FILE *log_file;
  struct object *child_objp;
  double tm[4][4];
  unsigned short l,m,n;
  char *tmp_name;
  int i;

  log_file=world->log_file;
  l=4;
  m=4;
  n=4;
  mult_matrix(objp->t_matrix,im,tm,l,m,n);

  if (sub_name!=NULL) { 
    if (strcmp(sub_name,"")==0) {
      tmp_name=my_strdup("");
      if(tmp_name == NULL){
		fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        	i = emergency_output();
		fprintf(stderr,"Fatal error: out of memory during bounding box computation.\nAttempt to write intermediate results had %d errors.\n", i);
        	exit(EXIT_FAILURE);
      }else{}
    }
    else {
      tmp_name=my_strcat(sub_name,".");              
      if(tmp_name == NULL){
		fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        	i = emergency_output();
		fprintf(stderr,"Fatal error: out of memory while bounding box computation.\nAttempt to write intermediate results had %d errors.\n", i);
        	exit(EXIT_FAILURE);
      }else{}
    }
    sub_name=my_strcat(tmp_name,objp->last_name);    
    if(sub_name == NULL){
		fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        	i = emergency_output();
		fprintf(stderr,"Fatal error: out of memory while bounding box computation.\nAttempt to write intermediate results had %d errors.\n", i);
        	exit(EXIT_FAILURE);
    }else{}
    free((void *)tmp_name);
  }
  else {
    sub_name=my_strdup(objp->last_name);    
    if(sub_name == NULL){
		fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        	i = emergency_output();
		fprintf(stderr,"Fatal error: out of memory while bounding box computation.\nAttempt to write intermediate results had %d errors.\n", i);
        	exit(EXIT_FAILURE);
    }else{}
  }

  switch (objp->object_type) {
  case META_OBJ:
    no_printf("Bounding box of Meta object %s is computed\n",sub_name);
    fflush(log_file);
    child_objp=objp->first_child;
    while (child_objp!=NULL) {
      if (compute_bb(child_objp,tm,sub_name)) {
        return(1);
      }
      child_objp=child_objp->next;
    }
    break;
  case REL_SITE_OBJ:
    no_printf("Bounding box of Release site object %s is computed\n",sub_name);
    fflush(log_file);
    if (compute_bb_release_site(objp,tm)) {
      return(1);
    }
    break;
  case BOX_OBJ:
    no_printf("Bounding box of Box object %s is computed\n",sub_name);
    fflush(log_file);
    if (compute_bb_polygon_object(objp,tm,sub_name)) {
      return(1);
    }
    break;
  case POLY_OBJ:
    no_printf("Bounding box of Polygon list object %s is computed\n",sub_name);
    fflush(log_file);
    if (compute_bb_polygon_object(objp,tm,sub_name)) {
      return(1);
    }
    break;
  }

  free((void *)sub_name);
  return(0);
}



/**
 * Updates the bounding box of the world based on the size
 * and location of a release site.
 * Used by compute_bb().
 */
int compute_bb_release_site(struct object *objp, double (*im)[4])
{
  struct release_site_obj *rsop;
  double location[1][4];
  unsigned short l,m,n;
  double diam_x, diam_y, diam_z; /* diameters of the release_site */ 
 
  rsop=(struct release_site_obj *)objp->contents;
  
  if (rsop->release_shape == SHAPE_REGION) return 0;

  l=1;
  m=4;
  n=4;
  location[0][0]=rsop->location->x;
  location[0][1]=rsop->location->y;
  location[0][2]=rsop->location->z;
  location[0][3]=1.0;
  mult_matrix(location,im,location,l,m,n);
  
  if(rsop->diameter == NULL){
	diam_x = diam_y = diam_z = 0;
  }else{
        diam_x = rsop->diameter->x;
        diam_y = rsop->diameter->y;
        diam_z = rsop->diameter->z;
  }


  if (location[0][0]  - diam_x < world->bb_min.x) {
    world->bb_min.x=location[0][0] - diam_x;
  }
  if (location[0][1] - diam_y < world->bb_min.y) {
    world->bb_min.y=location[0][1] - diam_y;
  }
  if (location[0][2] - diam_z < world->bb_min.z) {
    world->bb_min.z=location[0][2] - diam_z;
  }
  if (location[0][0] + diam_x > world->bb_max.x) {
    world->bb_max.x=location[0][0] + diam_x;
  }
  if (location[0][1] + diam_y > world->bb_max.y) {
    world->bb_max.y=location[0][1] + diam_y;
  }
  if (location[0][2] + diam_z > world->bb_max.z) {
    world->bb_max.z=location[0][2] + diam_z;
  }

  return(0);
}



/**
 * Updates the bounding box of the world based on the size
 * and location of a polygon_object.
 * Used by compute_bb().
 */
int compute_bb_polygon_object(struct object *objp, double (*im)[4], char *full_name)
{
  struct polygon_object *pop;
  struct ordered_poly *opp;
  double p[1][4];
  int i,n_verts,n_walls;
  unsigned short l,m,n;

  pop=(struct polygon_object *)objp->contents;
  n_walls=pop->n_walls;
  l=1;
  m=4;
  n=4;

    opp=(struct ordered_poly *)pop->polygon_data;
    n_verts=opp->n_verts;
 
    for (i=0;i<n_verts;i++) {
      p[0][0]=opp->vertex[i].x;
      p[0][1]=opp->vertex[i].y;
      p[0][2]=opp->vertex[i].z;
      p[0][3]=1.0;
      mult_matrix(p,im,p,l,m,n);
      if (p[0][0]/world->length_unit<world->bb_min.x) {
        world->bb_min.x=p[0][0]/world->length_unit;
      }
      if (p[0][1]/world->length_unit<world->bb_min.y) {
        world->bb_min.y=p[0][1]/world->length_unit;
      }
      if (p[0][2]/world->length_unit<world->bb_min.z) {
        world->bb_min.z=p[0][2]/world->length_unit;
      }
      if (p[0][0]/world->length_unit>world->bb_max.x) {
        world->bb_max.x=p[0][0]/world->length_unit;
      }
      if (p[0][1]/world->length_unit>world->bb_max.y) {
        world->bb_max.y=p[0][1]/world->length_unit;
      }
      if (p[0][2]/world->length_unit>world->bb_max.z) {
        world->bb_max.z=p[0][2]/world->length_unit;
      }
    }

/*
    fprintf(log_file,"MCell: Total area of physical object %s = %.9g microns^2\n",objp->sym->name,total_area*world->length_unit*world->length_unit);
*/

  return(0);
}



/**
 * Instantiates a polygon_object.
 * Creates walls from a template polygon_object or box object
 * as defined in the MDL file after applying the necessary geometric
 * transformations (scaling, rotation and translation).
 * <br>
 */
int instance_polygon_object(struct object *objp, double (*im)[4], struct viz_obj *vizp, struct lig_count_ref *obj_lcrp, char *full_name)
{
// #define INIT_VERTEX_NORMALS
// Uncomment to compute vertex normals
  FILE *log_file;
  struct polygon_object *pop;
  struct ordered_poly *opp;
  struct vector3 *v,**vp;
  struct wall *w,**wp;
  struct viz_child *vcp;
  double p[1][4];
#ifdef INIT_VERTEX_NORMALS
  struct vector3 *vertex_normal;
  double origin[1][4];
#endif
  double total_area;
  int i,n_verts,n_walls,index_0,index_1,index_2;
  unsigned int degenerate_count;
  unsigned short l,m,n;
  char *obj_name;
  byte compute_vertex_normals;

  log_file=world->log_file;
  pop=(struct polygon_object *)objp->contents;
  n_walls=pop->n_walls;
  n_verts=pop->n_verts;
  l=1;
  m=4;
  n=4;
  total_area=0;

  obj_name=my_strdup(full_name);
  if(obj_name == NULL)
  {
		fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        	int i = emergency_output();
		fprintf(stderr,"Fatal error: out of memory while instantion of polygon object.\nAttempt to write intermediate results had %d errors.\n", i);
        	exit(EXIT_FAILURE);
  }

/* Allocate and initialize walls and vertices */
    w=(struct wall *)malloc(n_walls*sizeof(struct wall));
    wp=(struct wall **)malloc(n_walls*sizeof(struct wall *));
    v=(struct vector3 *)malloc(n_verts*sizeof(struct vector3));
    vp=(struct vector3 **)malloc(n_verts*sizeof(struct vector3 *)); 
    if (w==NULL || wp==NULL || v==NULL || vp==NULL)
    {
      fprintf(world->err_file,"Out of memory while instantiating polygon object.  Quitting.\n");
      exit(EXIT_FAILURE);
    }
    objp->walls=w;
    objp->wall_p=wp;
    objp->verts=v;
    objp->vert_p=vp;

    opp=(struct ordered_poly *)pop->polygon_data;

    compute_vertex_normals=0;

/* If we want vertex normals we'll have to add a place to store them
   in struct object.
*/
#ifdef INIT_VERTEX_NORMALS
    if (opp->normal!=NULL) {
      compute_vertex_normals=1;
    }
#endif
   if(vizp!=NULL)
   {
     if((world->viz_mode == DREAMM_V3_MODE) || (world->viz_mode == DREAMM_V3_GROUPED_MODE) || (objp->viz_state!=NULL))
     {

      if ((vcp=(struct viz_child *)malloc
           (sizeof(struct viz_child)))==NULL) {
		fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        	int i = emergency_output();
		fprintf(stderr,"Fatal error: out of memory while instantiation of polygon object.\nAttempt to write intermediate results had %d errors.\n", i);
        	exit(EXIT_FAILURE);
      }
      vcp->obj = objp;
      vcp->next = vizp->viz_child_head;
      vizp->viz_child_head = vcp;
    }
  }  

  for (i=0;i<n_verts;i++) {
    vp[i]=&v[i];
    p[0][0]=opp->vertex[i].x;
    p[0][1]=opp->vertex[i].y;
    p[0][2]=opp->vertex[i].z;
    p[0][3]=1.0;
    mult_matrix(p,im,p,l,m,n);
    v[i].x=p[0][0]/world->length_unit;
    v[i].y=p[0][1]/world->length_unit;
    v[i].z=p[0][2]/world->length_unit;

#ifdef INIT_VERTEX_NORMALS
    if (compute_vertex_normals) {
      p[0][0]=opp->normal[i].x;
      p[0][1]=opp->normal[i].y;
      p[0][2]=opp->normal[i].z;
      p[0][3]=1.0;
      origin[0][0]=0;
      origin[0][1]=0;
      origin[0][2]=0;
      origin[0][3]=1.0;
      mult_matrix(p,im,p,l,m,n);
      mult_matrix(origin,im,origin,l,m,n);
      vertex_normal[i].x=p[0][0]-origin[0][0];
      vertex_normal[i].y=p[0][1]-origin[0][1];
      vertex_normal[i].z=p[0][2]-origin[0][2];
      normalize(&vertex_normal[i]);
    }
#endif
  }
  
  degenerate_count=0;
  for (i=0;i<n_walls;i++) {
    if (!get_bit(pop->side_removed,i)) {
      wp[i]=&w[i];
      index_0=opp->element[i].vertex_index[0];
      index_1=opp->element[i].vertex_index[1];
      index_2=opp->element[i].vertex_index[2];

      init_tri_wall(objp,i,vp[index_0],vp[index_1],vp[index_2]);
      total_area+=wp[i]->area;

      if (wp[i]->area==0) {
        if (world->notify->degenerate_polys != WARN_COPE)
        {
          if (world->notify->degenerate_polys==WARN_ERROR)
          {
            log_file = world->err_file;
            fprintf(log_file,"\nError -- ");
          }
          else fprintf(log_file,"\nWarning -- ");
          
          fprintf(log_file,"Degenerate polygon found and automatically removed: %s %d\n\n",objp->sym->name,i);
          fprintf(log_file,"  Vertex 0: %.5e %.5e %.5e\n",vp[index_0]->x,vp[index_0]->y,vp[index_0]->z);
          fprintf(log_file,"  Vertex 1: %.5e %.5e %.5e\n",vp[index_1]->x,vp[index_1]->y,vp[index_1]->z);
          fprintf(log_file,"  Vertex 2: %.5e %.5e %.5e\n",vp[index_2]->x,vp[index_2]->y,vp[index_2]->z);
          
          if (world->notify->degenerate_polys==WARN_ERROR) return 1;
        }
        set_bit(pop->side_removed,i,1);
        objp->n_walls_actual--;
        degenerate_count++;
        wp[i]=NULL;
      }
    }
    else {
      wp[i]=NULL;
    }
  }
  if (degenerate_count) remove_gaps_from_regions(objp);
  
  objp->total_area=total_area;
  
#ifdef DEBUG    
  printf("n_walls = %d\n", n_walls);
  printf("n_walls_actual = %d\n", objp->n_walls_actual);
#endif

  return(0);
}



int init_regions()
{
  if (world->clamp_list!=NULL) init_clamp_lists();

  if (instance_obj_regions(world->root_instance,NULL)) {
    return(1);
  }

  return(0);
}


/* First part of concentration clamp initialization. */
/* After this, list is grouped by surface class. */
/* Second part (list of objects) happens with regions. */
void init_clamp_lists()
{
  struct ccn_clamp_data *ccd,*temp;
  
  /* Sort by memory order of surface_class pointer--handy way to collect like classes */
  world->clamp_list = (struct ccn_clamp_data*)void_list_sort((struct void_list*)world->clamp_list);
  
  /* Toss other molecules in same surface class into next_mol lists */
  for (ccd = world->clamp_list ; ccd!=NULL ; ccd=ccd->next)
  {
    while (ccd->next != NULL && ccd->surf_class==ccd->next->surf_class)
    {
      ccd->next->next_mol = ccd->next_mol;
      ccd->next_mol = ccd->next;
      ccd->next = ccd->next->next;
    }
    for (temp=ccd->next_mol ; temp!=NULL ; temp=temp->next_mol)
    {
      temp->next = ccd->next;
    }
  }
}


int instance_obj_regions(struct object *objp,char *sub_name)
{
  FILE *log_file;
  struct object *child_objp;
  char *tmp_name;

  log_file=world->log_file;

  if (sub_name!=NULL) { 
    if (strcmp(sub_name,"")==0) {
      tmp_name=my_strdup("");
      if (tmp_name == NULL) {
		fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        	int i = emergency_output();
		fprintf(stderr,"Fatal error: out of memory while instantiation of object regions.\nAttempt to write intermediate results had %d errors.\n", i);
        	exit(EXIT_FAILURE);
      }else{}
    }
    else {
      tmp_name=my_strcat(sub_name,".");              
      if (tmp_name == NULL) {
		fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        	int i = emergency_output();
		fprintf(stderr,"Fatal error: out of memory while instantiation of object regions.\nAttempt to write intermediate results had %d errors.\n", i);
        	exit(EXIT_FAILURE);
      }else{}
    }
    sub_name=my_strcat(tmp_name,objp->last_name);    
    if (sub_name == NULL) {
		fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        	int i = emergency_output();
		fprintf(stderr,"Fatal error: out of memory while instantiation of object regions.\nAttempt to write intermediate results had %d errors.\n", i);
        	exit(EXIT_FAILURE);
    }else{}
    free((void *)tmp_name);
  }
  else {
    sub_name=my_strdup(objp->last_name);    
    if (sub_name == NULL) {
		fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        	int i = emergency_output();
		fprintf(stderr,"Fatal error: out of memory while instantiation of object regions.\nAttempt to write intermediate results had %d errors.\n", i);
        	exit(EXIT_FAILURE);
    }else{}
  }

  switch (objp->object_type) {
  case META_OBJ:
    no_printf("Initializing regions in meta object: %s\n",sub_name);
    fflush(log_file);
    child_objp=objp->first_child;
    while (child_objp!=NULL) {
      if (instance_obj_regions(child_objp,sub_name)) {
        return(1);
      }
      child_objp=child_objp->next;
    }
    break;
  case REL_SITE_OBJ:
    break;
  case BOX_OBJ:
    no_printf("Initializing regions in box object: %s\n",sub_name);
    fflush(log_file);
    if (init_wall_regions(objp,sub_name)) {
      return(1);
    }
    break;
  case POLY_OBJ:
    no_printf("Initializing regions in polygon list object: %s\n",sub_name);
    fflush(log_file);
    if (init_wall_regions(objp,sub_name)) {
      return(1);
    }
    break;
  }

  free((void *)sub_name);
  return(0);
}



/**
 * Initialize data associated with wall regions.
 * This function is called during wall instantiation Pass #3
 * after walls have been copied to sub-volume local memory.
 * Sets wall surf_class by region.
 * Creates surface grids.
 * Populates effector tiles by region.
 * Creates virtual regions on which to clamp concentration
 */
int init_wall_regions(struct object *objp, char *full_name)
{
  FILE *log_file;
  struct polygon_object *pop;
  struct wall *w;
  struct region *rp;
  struct region_list *rlp,*wrlp;
  int i,n_walls;

  log_file=world->log_file;

  pop=(struct polygon_object *)objp->contents;
  n_walls=pop->n_walls;

   
  no_printf("Processing %d regions in polygon list object: %s\n",objp->num_regions,full_name);

  /* prepend a copy of eff_dat for each element referenced in each region
     of this object to the eff_prop list for the referenced element */
  rlp=objp->regions;
  
  for (rlp=objp->regions ; rlp!=NULL ; rlp=rlp->next)
  {
    rp=rlp->reg;
    if (rp->membership==NULL)
    {
      fprintf(world->err_file,"Internal error: incomplete region information for %s\n",rp->sym->name);
      return 1;
    }
    for (i=0;i<rp->membership->nbits;i++)
    {
      if (get_bit(rp->membership,i))
      {
	/* prepend this region to wall region list of i_th wall only if the region is used in counting */
	w=objp->wall_p[i];
	rp->area += w->area;
	if (rp->surf_class!=NULL) w->surf_class = rp->surf_class;

	if ((rp->flags & COUNT_SOME) != 0)
	{  
	  wrlp = (struct region_list *)mem_get(w->birthplace->regl);
	  if (wrlp==NULL)
	  {
	    fprintf(world->err_file,"Out of memory: can't place regions on geometry.\n");
	    return 1;
	  }
	  wrlp->reg=rp;
	  wrlp->next=w->counting_regions;
	  w->counting_regions=wrlp;
	  w->flags|=rp->flags;
	}
      }
    }
  } /*end loop over all regions in object */
  no_printf("Total area of object %s = %.9g um^2\n",objp->sym->name,objp->total_area/world->effector_grid_density);
  no_printf("  number of tiles = %u\n",objp->n_tiles);
  no_printf("  number of occupied tiles = %u\n",objp->n_occupied_tiles);
  no_printf("  grid molecule density = %.9g\n",objp->n_occupied_tiles*world->effector_grid_density/objp->total_area);
    
  /* Check to see if we need to generate virtual regions for */
  /* concentration clamps on this object */
  if (world->clamp_list!=NULL)
  {
    struct ccn_clamp_data *ccd;
    struct ccn_clamp_data *temp;
    int j;
    int found_something = 0;
    
    for (i=0;i<n_walls;i++)
    {
      if (get_bit(pop->side_removed,i)) continue;
      if (objp->wall_p[i]->surf_class != world->g_surf)
      {
        for (ccd=world->clamp_list ; ccd!=NULL ; ccd=ccd->next)
        {
          if (objp->wall_p[i]->surf_class == ccd->surf_class)
          {
            if (ccd->objp!=objp)
            {
              if (ccd->objp==NULL) ccd->objp=objp;
              else if (ccd->next_obj != NULL && ccd->next_obj->objp==objp) ccd=ccd->next_obj;
              else
              {
                temp = (struct ccn_clamp_data*)malloc(sizeof(struct ccn_clamp_data));
                if (temp==NULL)
                {
                  fprintf(world->err_file,"Out of memory assembling concentration clamp data.\n");
                  return 1;
                }
                memcpy(temp,ccd,sizeof(struct ccn_clamp_data));
                temp->objp = objp;
                ccd->sides = NULL;
                ccd->n_sides = 0;
                ccd->side_idx = NULL;
                ccd->cum_area = NULL;
                ccd->next_obj = temp;
                ccd = temp;
              }
            }
            if (ccd->sides==NULL)
            {
              ccd->sides = new_bit_array(n_walls);
              if (ccd->sides==NULL)
              {
                fprintf(world->err_file,"Out of memory assembling concentration clamp data.\n");
                return 1;
              }
              set_all_bits(ccd->sides,0);
            }
            set_bit(ccd->sides,i,1);
            ccd->n_sides++;
            found_something=1;
          }
        }
      }
    }
    
    if (found_something)
    {
      for (ccd=world->clamp_list ; ccd!=NULL ; ccd=ccd->next)
      {
        if (ccd->objp!=objp)
        {
          if (ccd->next_obj!=NULL && ccd->next_obj->objp==objp) ccd=ccd->next_obj;
          else continue;
        }
        
        ccd->side_idx = (int*)malloc(ccd->n_sides*sizeof(int));
        ccd->cum_area = (double*)malloc(ccd->n_sides*sizeof(double));
        if (ccd->side_idx==NULL || ccd->cum_area==NULL)
        {
          fprintf(world->err_file,"Out of memory assembling concentration clamp data.\n");
          return 1;
        }
        
        j=0;
        for (i=0;i<n_walls;i++)
        {
          if (get_bit(ccd->sides,i))
          {
            ccd->side_idx[j] = i;
            ccd->cum_area[j] = objp->wall_p[i]->area;
            j++;
          }
        }
        if (j!=ccd->n_sides)
        {
          fprintf(world->err_file,"Miscounted the number of walls for concentration clamp\n  on object %s\n  surface class %s\n",
                  objp->sym->name,ccd->surf_class->sym->name);
          return 1;
        }
        
        for (j=1;j<ccd->n_sides;j++) ccd->cum_area[j] += ccd->cum_area[j-1];
        
        ccd->scaling_factor = ccd->cum_area[ccd->n_sides-1] * world->length_unit * world->length_unit * world->length_unit /
                              2.9432976599069717358e-9;  /* sqrt(MY_PI)/(1e-15*N_AV) */
        if (ccd->orient!=0) ccd->scaling_factor *= 0.5;
      }
    }
    
  }


#ifdef KELP
  cdp->sym->ref_count--;
  if (!cdp->sym->ref_count) {	/* Done with the geometry information */
	destroy_sym_value(cdp->sym);	/* free up memory */
  }
#endif

  return(0);
}



int init_effectors()
{
  if (instance_obj_effectors(world->root_instance)) return 1;
  return 0;
}


int instance_obj_effectors(struct object *objp)
{
  struct object *child_objp;

  switch (objp->object_type)
  {
    case META_OBJ:
      for (child_objp=objp->first_child ; child_objp!=NULL ; child_objp=child_objp->next)
      {
	if (instance_obj_effectors(child_objp)) return 1;
      }
      break;
    case REL_SITE_OBJ:
      break;
    case BOX_OBJ:
    case POLY_OBJ:
      if (init_wall_effectors(objp)) return 1;
      break;
    default:
      break;
  }

  return(0);
}


int init_wall_effectors(struct object *objp)
{
  FILE *log_file;
  struct polygon_object *pop;
  struct wall *w;
  struct eff_dat *effdp,*dup_effdp,**eff_prop;
  struct region *rp;
  struct region_list *rlp,*rlp2,*reg_eff_num_head;
  int i,n_walls;
  byte reg_eff_num;

  log_file=world->log_file;

  pop=(struct polygon_object *)objp->contents;
  n_walls=pop->n_walls;

   
  /* allocate scratch storage to hold effector info for each wall */
  if ((eff_prop=(struct eff_dat **)malloc(n_walls*sizeof(struct eff_dat *)))==NULL)
  {
    fprintf(world->err_file,"Out of memory: can't create space for molecules on a region.\n");
    return 1;
  }

  for (i=0;i<n_walls;i++) eff_prop[i]=NULL; 

  /* prepend a copy of eff_dat for each element referenced in each region
     of this object to the eff_prop list for the referenced element */
  reg_eff_num_head=NULL;
  rlp=objp->regions;
  
  for (rlp=objp->regions ; rlp!=NULL ; rlp=rlp->next)
  {
    rp=rlp->reg;
    reg_eff_num=0;

    for (i=0;i<rp->membership->nbits;i++)
    {
      if (get_bit(rp->membership,i))
      {
	w=objp->wall_p[i];

	/* prepend region eff data for this region to eff_prop for i_th wall */
        for ( effdp=rp->eff_dat_head ; effdp!=NULL ; effdp=effdp->next )
	{
	  if (effdp->quantity_type==EFFDENS)
	  {
	    if ((dup_effdp=(struct eff_dat *)malloc(sizeof(struct eff_dat)))==NULL)
	    {
	      fprintf(world->err_file,"Out of memory: can't create space for molecules on a region.\n");
	      return 1;
	    }
	    dup_effdp->eff=effdp->eff;
	    dup_effdp->quantity_type=effdp->quantity_type;
	    dup_effdp->quantity=effdp->quantity;
	    dup_effdp->orientation=effdp->orientation;
	    dup_effdp->next=eff_prop[i];
	    eff_prop[i]=dup_effdp;
	  }
	  else reg_eff_num=1;
	}

	/* prepend surf_class eff data for this region to eff_prop for i_th wall on last region */
	if (w->surf_class != NULL && rlp->next==NULL)
	{
	  for ( effdp=w->surf_class->eff_dat_head ; effdp!=NULL ; effdp=effdp->next )
	  {
	    if (effdp->quantity_type==EFFDENS)
	    {
	      if ((dup_effdp=(struct eff_dat *)malloc(sizeof(struct eff_dat)))==NULL)
	      {
		fprintf(world->err_file,"Out of memory: can't create space for molecules on a region.\n");
		return 1;
	      }
	      dup_effdp->eff=effdp->eff;
	      dup_effdp->quantity_type=effdp->quantity_type;
	      dup_effdp->quantity=effdp->quantity;
	      dup_effdp->orientation=effdp->orientation;
	      dup_effdp->next=eff_prop[i];
	      eff_prop[i]=dup_effdp;
	    }
	    else reg_eff_num=1;
	  }
	}

      }
    } /* done checking each wall */
    
    if (reg_eff_num)
    {
      if ((rlp2=(struct region_list *)malloc(sizeof(struct region_list)))==NULL)
      {
	fprintf(world->err_file,"Out of memory: can't place regions on geometry.\n");
	return 1;
      }
      rlp2->reg=rp;
      rlp2->next=reg_eff_num_head;
      reg_eff_num_head=rlp2;
    }
  } /*end for (... ; rlp != NULL ; ...) */

  for (i=0;i<n_walls;i++)
  {
    if (!get_bit(pop->side_removed,i))
    {
      if (eff_prop[i]!=NULL)
      {
	if (init_effectors_by_density(objp->wall_p[i],eff_prop[i])) return 1;
      }
    }
  }

  if (reg_eff_num_head!=NULL)
  {
    if (init_effectors_by_number(objp,reg_eff_num_head)) {
      return(1);
    }
    
    /* free region list created to hold regions populated by number */
    rlp=reg_eff_num_head;
    while(rlp!=NULL)
    {
      rlp2=rlp;
      rlp=rlp->next;
      free(rlp2);
    }
  }

  /* free eff_prop array and contents */
  for (i=0;i<n_walls;i++)
  {
    if (eff_prop[i]!=NULL)
    {
      effdp=eff_prop[i];
      while(effdp!=NULL)
      {
	dup_effdp=effdp;
	effdp=effdp->next;
	free(dup_effdp);
      }
    }
  }
  free(eff_prop);

    
  return(0);
}


int init_effectors_by_density(struct wall *w, struct eff_dat *effdp_head)
{
  FILE *log_file;
  struct object *objp;
  struct species **eff;
  struct surface_grid *sg;
  struct eff_dat *effdp;
  struct grid_molecule *mol;
  short *orientation;
  unsigned int i,j,n,nr,n_occupied;
  int p_index;
  double rand,*prob,area,tot_prob,tot_density;

  log_file=world->log_file;

  no_printf("Initializing effectors by density...\n");
  fflush(log_file);

  if (create_grid(w,NULL)) {
    return(1);
  }
  sg=w->effectors;
  objp=w->parent_object;

  nr=0;
  effdp=effdp_head;
  while (effdp!=NULL) {
    nr++;
    effdp=effdp->next;
  }

  if ((eff=(struct species **)malloc(nr*sizeof(struct species *)))==NULL) {
			fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        		int i = emergency_output();
			fprintf(stderr,"Fatal error: out of memory while effectors by density intialization.\nAttempt to write intermediate results had %d errors.\n", i);
        		exit(EXIT_FAILURE);
  }
  if ((prob=(double *)malloc(nr*sizeof(double)))==NULL) {
			fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        		int i = emergency_output();
			fprintf(stderr,"Fatal error: out of memory while effectors by density intialization.\nAttempt to write intermediate results had %d errors.\n", i);
        		exit(EXIT_FAILURE);
  }
  if ((orientation=(short*)malloc(nr*sizeof(short)))==NULL) {
			fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        		int i = emergency_output();
			fprintf(stderr,"Fatal error: out of memory while effectors by density intialization.\nAttempt to write intermediate results had %d errors.\n", i);
        		exit(EXIT_FAILURE);
  }

  for (i=0;i<nr;i++) {
    eff[i]=NULL;
    prob[i]=0.0;
    orientation[i]=0;
  }

  n=sg->n_tiles;
  area=w->area;
  objp->n_tiles+=n;
  no_printf("Initializing %d effectors...\n",n);
  no_printf("  Area = %.9g\n",area);
  no_printf("  Grid_size = %d\n",sg->n);
  no_printf("  Number of effector types in wall = %d\n",nr);

  i=0;
  tot_prob=0;
  tot_density=0;
  effdp=effdp_head;
  while (effdp!=NULL) {
    no_printf("  Adding effector %s to wall at density %.9g\n",effdp->eff->sym->name,effdp->quantity);
    tot_prob+=(area*effdp->quantity)/(n*world->effector_grid_density);
    prob[i]=tot_prob;
    if (effdp->orientation > 0) orientation[i] = 1;
    else if (effdp->orientation < 0) orientation[i] = -1;
    else orientation[i] = (rng_uint(world->rng)&1)?1:-1;
    eff[i++]=effdp->eff;
    tot_density+=effdp->quantity;
    effdp=effdp->next;
  }

  if (tot_density>world->effector_grid_density) {
    fprintf(log_file,"\nMCell: Warning -- Total effector density too high: %f\n\n",tot_density);
    fflush(log_file);
/*
    return(1);
*/
  }

  n_occupied=0;
  for (i=0;i<n;i++) {
    if (world->chkpt_init) {
      j=0;
      p_index=-1;
      rand = rng_dbl(world->rng);
      while (j<nr && p_index==-1) {
        if (rand<=prob[j++]) {
          p_index=j-1;
        }
      }
      if (p_index!=-1) {
        n_occupied++;
        eff[p_index]->population++;
        mol=(struct grid_molecule *)mem_get(w->birthplace->gmol);
        if(mol == NULL){
			fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        		int i = emergency_output();
			fprintf(stderr,"Fatal error: out of memory while effectors by density intialization.\nAttempt to write intermediate results had %d errors.\n", i);
        		exit(EXIT_FAILURE);
        }
        sg->mol[i]=mol;
        mol->t=0;
        mol->t2=0;
	mol->birthday = 0;

        mol->properties=eff[p_index];
        mol->birthplace=w->birthplace->gmol;
        mol->grid_index=i;
	if (world->randomize_gmol_pos) grid2uv_random(sg,i,&(mol->s_pos));
	else grid2uv(sg,i,&(mol->s_pos));
        mol->orient=orientation[p_index];
        mol->grid=sg;

        mol->flags=TYPE_GRID|ACT_NEWBIE|IN_SCHEDULE|IN_SURFACE;
	if (mol->properties->space_step>0) mol->flags |= ACT_DIFFUSE;
        if ( trigger_unimolecular(eff[p_index]->hashval,(struct abstract_molecule *)mol)!=NULL
	     || (eff[p_index]->flags&CAN_GRIDWALL)!=0 ) {
          mol->flags|=ACT_REACT;
        }

        if ((mol->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED)) != 0)
          count_me_by_region((struct abstract_molecule*)mol,1,NULL);
      
        if (schedule_add(w->birthplace->timer,mol)){ 
		fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        	int i = emergency_output();
		fprintf(stderr,"Fatal error: out of memory while effectors by density intialization.\nAttempt to write intermediate results had %d errors.\n", i);
        	exit(EXIT_FAILURE);
         }
      }
    }
  }

  sg->n_occupied=n_occupied;
  objp->n_occupied_tiles+=n_occupied;

  for (i=0;i<nr;i++) {
    no_printf("Total number of effector %s = %d\n",eff[i]->sym->name,eff[i]->population);
  }

  free(eff);
  free(prob);
  free(orientation);

  no_printf("Done initializing %d effectors by density\n",n_occupied);
  fflush(log_file);

	
  return(0);
}



int init_effectors_by_number(struct object *objp, struct region_list *reg_eff_num_head)
{
  FILE *log_file;
  struct polygon_object *pop;
  struct species *eff;
  struct grid_molecule ***tiles,***tiles_tmp;
  struct grid_molecule gmol,*bread_crumb,*mol;
  struct region_list *rlp; 
  struct region *rp;
/*  struct element_list *elp; */
  struct surface_grid *sg;
  struct eff_dat *effdp;
  struct wall **walls,**walls_tmp,*w;
  short orientation;
  unsigned int *index,*index_tmp;
  unsigned int n_free_eff,n_set,n_clear;
  unsigned int i,j,k;
  byte done;

    no_printf("Initializing effectors by number...\n");

    log_file=world->log_file;
    pop=(struct polygon_object *)objp->contents;

    tiles=NULL;
    tiles_tmp=NULL;
    index=NULL;
    index_tmp=NULL;
    walls=NULL;
    walls_tmp=NULL;
    bread_crumb=&gmol;

    /* traverse region list and add effector sites by number to whole regions
       as appropriate */
    rlp=reg_eff_num_head;
    while (rlp!=NULL) {
      rp=rlp->reg;
        /* initialize effector grids in region as needed and */
        /* count total number of free effector sites in region */
        n_free_eff=0;
	for (i=0;i<rp->membership->nbits;i++)
	{
	  if (get_bit(rp->membership,i))
	  {
	    w=objp->wall_p[i];
	    if (create_grid(w,NULL)) {
	      return(1);
	    }
	    sg=w->effectors;
	    n_free_eff=n_free_eff+(sg->n_tiles-sg->n_occupied);
          }
        }
        no_printf("Number of free effector tiles in region %s = %d\n",rp->sym->name,n_free_eff);
        fflush(stdout);
      if (world->chkpt_init) {  /* only needed for denovo initiliazation */
        /* allocate memory to hold array of pointers to all free tiles */
        if ((tiles=(struct grid_molecule ***)malloc
           (n_free_eff*sizeof(struct grid_molecule **)))==NULL) {
			fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        		int i = emergency_output();
			fprintf(stderr,"Fatal error: out of memory while effectors by number intialization.\nAttempt to write intermediate results had %d errors.\n", i);
        		exit(EXIT_FAILURE);
        }
        if ((index=(unsigned int *)malloc
           (n_free_eff*sizeof(unsigned int)))==NULL) {
			fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        		int i = emergency_output();
			fprintf(stderr,"Fatal error: out of memory while effectors by number intialization.\nAttempt to write intermediate results had %d errors.\n", i);
        		exit(EXIT_FAILURE);
        }
        if ((walls=(struct wall **)malloc
           (n_free_eff*sizeof(struct wall *)))==NULL) {
			fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        		int i = emergency_output();
			fprintf(stderr,"Fatal error: out of memory while effectors by number intialization.\nAttempt to write intermediate results had %d errors.\n", i);
        		exit(EXIT_FAILURE);
        }
        /* initialize array of pointers to all free tiles */
        k=0;
	for (i=0;i<rp->membership->nbits;i++)
	{
	  if (get_bit(rp->membership,i))
	  {
	    w=objp->wall_p[i];
	    sg=w->effectors;
	    if (sg!=NULL) {
	      for (j=0;j<sg->n_tiles;j++) {
		if (sg->mol[j]==NULL) {
		  tiles[k]=&(sg->mol[j]);
		  index[k]=j;
		  walls[k++]=w;
		}
	      }
	    }
	  }
        }
      } /* end while(world->chkpt_init) */

      /* distribute desired number of effector sites */
      /* for each effector type to add */
      effdp=rp->eff_dat_head;
      while (effdp!=NULL) {
        if (effdp->quantity_type==EFFNUM) {
          eff=effdp->eff;

          if (world->chkpt_init) {  /* only needed for denovo initiliazation */
	    if (effdp->orientation > 0) orientation = 1;
	    else if (effdp->orientation < 0) orientation = -1;
	    else orientation = (rng_uint(world->rng)&1)?1:-1;
  
            n_set=effdp->quantity;
            n_clear=n_free_eff-n_set;
            eff->population+=n_set;

            if (n_set > n_free_eff) {
              fprintf(log_file,"\nMCell: Warning -- Number of %s effectors to place (%d) exceeds number of free effector tiles (%d) in region %s[%s]\n\n",eff->sym->name,n_set,n_free_eff,rp->parent->sym->name,rp->region_last_name);
              n_set=n_free_eff;
              n_clear=0;
            }
            no_printf("distribute %d of effector %s\n",n_set,eff->sym->name);
            no_printf("n_set = %d  n_clear = %d  n_free_eff = %d\n",n_set,n_clear,n_free_eff);
            fflush(stdout);

            /* if filling more than half the free tiles
              init all with bread_crumbs
              choose which tiles to free again
              and then convert remaining bread_crumbs to actual molecules */
            if (n_set > n_free_eff/2) {
              no_printf("filling more than half the free tiles: init all with bread_crumb\n");
              fflush(stdout);
              for (j=0;j<n_free_eff;j++) {
                *tiles[j]=bread_crumb;
              }

              no_printf("choose which tiles to free again\n");
              fflush(stdout);
              for (j=0;j<n_clear;j++) {
                done=0;
                while (!done) {
		  k = (int) (rng_dbl(world->rng)*n_free_eff);
                  if (*tiles[k]==bread_crumb) {
                    *tiles[k]=NULL;
                    done=1;
                  }
                }
              }

              no_printf("convert remaining bread_crumbs to actual molecules\n");
              fflush(stdout);
              for (j=0;j<n_free_eff;j++) {
                if (*tiles[j]==bread_crumb) {
                  mol=(struct grid_molecule *)
                    mem_get(walls[j]->birthplace->gmol);
                  if (mol == NULL){
			fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        		int i = emergency_output();
			fprintf(stderr,"Fatal error: out of memory while effectors by number initialization.\nAttempt to write intermediate results had %d errors.\n", i);
        		exit(EXIT_FAILURE);
                  }
                  *tiles[j]=mol;
                  mol->t=0;
                  mol->t2=0;
		  mol->birthday=0;
                  mol->properties=eff;
                  mol->birthplace=walls[j]->birthplace->gmol;
                  mol->grid_index=index[j];
		  if (world->randomize_gmol_pos) grid2uv_random(walls[j]->effectors,index[j],&(mol->s_pos));
		  else grid2uv(walls[j]->effectors,index[j],&(mol->s_pos));
                  mol->orient=orientation;
                  mol->grid=walls[j]->effectors;
                  mol->flags=TYPE_GRID|ACT_NEWBIE|IN_SCHEDULE|IN_SURFACE;
		  if (mol->properties->space_step > 0) mol->flags |= ACT_DIFFUSE;
                  if (trigger_unimolecular(eff->hashval,(struct abstract_molecule *)mol)!=NULL
		      || (eff->flags&CAN_GRIDWALL)!=0 ) {
                    mol->flags|=ACT_REACT;
                  }
                  
                  if ((mol->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED)) != 0)
                    count_me_by_region((struct abstract_molecule*)mol,1,NULL);
      
                  if ( schedule_add(walls[j]->birthplace->timer,mol) ){ 
			fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        		int i = emergency_output();
			fprintf(stderr,"Fatal error: out of memory while effectors by number initialization.\nAttempt to write intermediate results had %d errors.\n", i);
        		exit(EXIT_FAILURE);
		   }
                }
              }
            }
            else {  /* just fill only the tiles we need */
              no_printf("fill only the tiles we need\n");
              fflush(stdout);
              for (j=0;j<n_set;j++) {
                done=0;
                while (!done) {
		  k = (int) (rng_dbl(world->rng)*n_free_eff);
                  if (*tiles[k]==NULL) {
                    mol=(struct grid_molecule *)mem_get
                      (walls[k]->birthplace->gmol);
                    if (mol == NULL){
			fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        		int i = emergency_output();
			fprintf(stderr,"Fatal error: out of memory while effectors by number initialization.\nAttempt to write intermediate results had %d errors.\n", i);
        		exit(EXIT_FAILURE);
                    }
                    *tiles[k]=mol;
                    mol->t=0;
                    mol->t2=0;
		    mol->birthday=0;
                    mol->properties=eff;
                    mol->birthplace=walls[k]->birthplace->gmol;
                    mol->grid_index=index[k];
		    if (world->randomize_gmol_pos) grid2uv_random(walls[k]->effectors,index[k],&(mol->s_pos));
		    else grid2uv(walls[k]->effectors,index[k],&(mol->s_pos));
                    mol->orient=orientation;
                    mol->grid=walls[k]->effectors;
                    mol->flags=TYPE_GRID|ACT_NEWBIE|IN_SCHEDULE|IN_SURFACE;
		    if (mol->properties->space_step > 0) mol->flags |= ACT_DIFFUSE;
                      if (trigger_unimolecular(eff->hashval,(struct abstract_molecule *)mol)!=NULL
		          || (eff->flags&CAN_GRIDWALL)!=0) {
                      mol->flags|=ACT_REACT;
                    }
                  
                    if ((mol->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED)) != 0)
                      count_me_by_region((struct abstract_molecule*)mol,1,NULL);
      
                    if ( schedule_add(walls[k]->birthplace->timer,mol) ){ 
			fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        		int i = emergency_output();
			fprintf(stderr,"Fatal error: out of memory while effectors by number initialization.\nAttempt to write intermediate results had %d errors.\n", i);
        		exit(EXIT_FAILURE);
                     }
                    done=1;
                  }
                }
              }
            }
        /* allocate memory to hold array of pointers to remaining free tiles */
            if ((tiles_tmp=(struct grid_molecule ***)malloc
                 (n_clear*sizeof(struct grid_molecule **)))==NULL) {
			fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        		int i = emergency_output();
			fprintf(stderr,"Fatal error: out of memory while effectors by number initialization.\nAttempt to write intermediate results had %d errors.\n", i);
        		exit(EXIT_FAILURE);
            }
            if ((index_tmp=(unsigned int *)malloc
               (n_clear*sizeof(unsigned int)))==NULL) {
			fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        		int i = emergency_output();
			fprintf(stderr,"Fatal error: out of memory while effectors by number initialization.\nAttempt to write intermediate results had %d errors.\n", i);
        		exit(EXIT_FAILURE);
            }
            if ((walls_tmp=(struct wall **)malloc
               (n_clear*sizeof(struct wall *)))==NULL) {
			fprintf(stderr,"Out of memory:trying to save intermediate results.\n");
        		int i = emergency_output();
			fprintf(stderr,"Fatal error: out of memory while effectors by number initialization.\nAttempt to write intermediate results had %d errors.\n", i);
        		exit(EXIT_FAILURE);
            }
            k=0;
            for (i=0;i<n_free_eff;i++) {
              if (*tiles[i]==NULL) {
                tiles_tmp[k]=tiles[i];
                index_tmp[k]=index[i];
                walls_tmp[k++]=walls[i];
              }
            }
            /* free original array of pointers to all free tiles */
            free(tiles);
            free(index);
            free(walls);
            tiles=tiles_tmp;
            index=index_tmp;
            walls=walls_tmp;
            n_free_eff=n_free_eff-n_set;

            /* update n_occupied for each effector grid */
	    for (i=0;i<rp->membership->nbits;i++)
	    {
	      if (get_bit(rp->membership,i))
	      {
		sg=objp->wall_p[i]->effectors;
		if (sg!=NULL)
		{
		  sg->n_occupied=0;
		  for (j=0;j<sg->n_tiles;j++) {
		    if (sg->mol[j]!=NULL) {
		      sg->n_occupied++;
		    }
		  }
                }
              }
            }
          } /* end while(world->chkpt_init) */
        }
        effdp=effdp->next;
      }
      /* free array of pointers to all free tiles */
      if (tiles!=NULL) {
        free(tiles);
      }
      if (index!=NULL) {
        free(index);
      }
      if (walls!=NULL) {
        free(walls);
      }
      rlp=rlp->next;
    }
    no_printf("Done initialize effectors by number.\n");
    return(0);
}






/* **************************************************************** */

# if 0






/**
 ** Initialize region counters for the RX_STATE and save 
 ** reg_counter_ref_list to the counter_hash_table of the object
 ** Start from MCell 2.68
 */
int init_region_counter(struct polygon_object *pop, struct cmprt_data *cdp, struct region_list *reg_count_head)
{
  struct rx *rx;
  struct region_list *rlp; 
  struct region *rp;
  struct element_list *elp;
  struct effector *ep;
  struct reg_counter_ref_list *rcrlp;
  struct reg_counter_ref *rcrp;
  int i,j;
  
    /* traverse region list and add effector sites by number to whole regions
       as appropriate */
    rlp=reg_count_head;
    no_printf(" Initialzing region counter ...\n");
    fflush(log_file);
    while (rlp!=NULL) {
      rp=rlp->region;
      elp=rp->element_list_head;
      /* check effector numbers for each region of the object*/
      while (elp!=NULL) {
	for (i=elp->begin;i<=elp->end;i++) {
	  if (pop->side_stat[i]) {
	    ep=cdp->wall[pop->cmprt_side_map[i]]->effectors;
	    if (ep!=NULL) {
	      for (j=0;j<ep->n_tiles;j++) {
		rx=ep->tiles[j];
		rcrlp=rp->reg_counter_ref_list;
		while (rcrlp!=NULL) {
		  rcrp=rcrlp->reg_counter_ref;
		  if ((rx!=NULL)&&(rx==rcrp->state)&&(rcrp->count_type==RX_STATE)) {
		    rcrp->counter++;
		  }
		  rcrlp=rcrlp->next;
		}
	      }
	    }
	  }
	}
	elp=elp->next;
      }
      rlp=rlp->next;
    }
    no_printf(" Done initializing region counter\n");
    fflush(log_file);
  return(0);
}


/* Initialize the region_list that contains this wall 
 * Add from MCell 2.68, designed for the region counters.
 */
struct region_list  *init_region_list_for_wall(struct region_list *rlp, int index) {
/*  struct element_list *elp; */
  struct region *rp;
  struct region_list *rlp1, *wall_region_head;
  
  wall_region_head=NULL;
  
  for ( ; rlp!=NULL ; rlp=rlp->next )
  {
    if (get_bit(rlp->region->membership,index)) {
      if ((rlp1=(struct region_list *)malloc(sizeof(struct region_list)))==NULL) {
	mdlerror("Can not save region_list for wall");
	return(NULL);
      }
      rlp1->region=rp;
      rlp1->next=wall_region_head;
      wall_region_head=rlp1;
    }
  }
  return(wall_region_head);

  /* free memory not needed any more*/  
  if (wall_region_head!=NULL) {
    rlp=wall_region_head;
    while (rlp!=NULL) {
      rlp1=rlp;
      rlp=rlp->next;
      free(rlp1);
    }
  }
}


/*
   ** Initialize the counter_hash_table of the object if 
   ** there are any reg_counters defined in this object, 
   ** and save the reg_counter_ref_list to the table.
   ** The table will be used for future counter searching.
   ** Start from MCell 2.68
 */  

int init_counter_hash_table (struct object *objp, struct region_list *reg_count_head) {

  struct counter_hash_table **countertab;
  struct counter_hash_table *chtp,*prev,*hash_table_head;
  struct region_list *rlp;
  struct region *rp;
  struct reg_counter_ref_list *rcrlp, *rcrlp2, *rcrlp_tmp;
  struct reg_counter_ref *rcrp, *rcrp2;
  char *rx_name;
  int col, done, i;
  unsigned short hashval;

  rlp=reg_count_head;
  rcrlp2=NULL;
  rcrp2=NULL;
  hash_table_head=NULL;
  no_printf("\t Initializing counter hash table for object %s\n",objp->sym->name);
  fflush(log_file);
  /* Initialize the object counter hash table if it is empty*/
  if (objp->counter_hash_table==NULL) {
    if ((countertab=(struct counter_hash_table **)malloc(COUNTER_HASH*sizeof(struct counter_hash_table *)))==NULL) {
      mdlerror("Can not save counter table to object");
      return(1);
    }
    countertab=init_countertab(COUNTER_HASH);
    objp->counter_hash_table=countertab;
  }
  
  /* save each reg_counter_ref_list to the table by the hash value of 
   * the rx state of the counter and sort by the address of the region
   */  
  while (rlp!=NULL) {
    rp=rlp->region;
    rcrlp=rp->reg_counter_ref_list;
    while (rcrlp!=NULL) {
      rcrp=rcrlp->reg_counter_ref;
      rx_name=rcrp->state->sym->name;
      hashval=hash(rx_name)&0x0000000f;
      /* Store counter in table, but check if it is already saved */
      if ((chtp=retrieve_counter(rx_name,rcrp, objp->counter_hash_table))==NULL) {
       if ((chtp=store_counter(rx_name,rcrlp, objp->counter_hash_table))==NULL) {
	fprintf(log_file,"Cannot store counter in table: %s\n",rx_name);
	return(NULL);
      }
      }
      else {
	warning(" Counter already defined");
	return(NULL);
      }
      objp->counter_hash_table[hashval]=chtp;      
      rcrlp=rcrlp->next;
    }
  rlp=rlp->next;
  }
  /*Print the counter_hash_table for debuging*/
  /*  
  for (i=0; i<COUNTER_HASH; i++) {
    hash_table_head=objp->counter_hash_table[i];
    while (hash_table_head!=NULL) {
      rcrlp2=(struct reg_counter_ref_list *)hash_table_head->value;
      if (rcrlp2!=NULL) {
	rcrp2=rcrlp2->reg_counter_ref;
	no_printf("\t\t col %d, counter %s, region %s = %d, type %d, counter %d;\n",i,rcrp2->state->sym->name,rcrp2->parent->sym->name,(int)rcrp2->parent, rcrp2->count_type,rcrp2->counter);
	fflush(log_file);
       }
      hash_table_head=hash_table_head->next;
     }
  }
  */  	
  no_printf("\t end of counter hash table initializing\n");
  fflush(log_file); 
  return(0);
}



int init_effector_grid(struct wall *wp)
{
  struct effector *ep;
  struct rx **tiles;
  struct vector3 v1,u_axis,v_axis,p_b,p_c,p_d,i_axis,j_axis;
  struct vector3 ab,bc,ac,step_u,step_v,diagonal,p0,p1,p2;
  unsigned short *dsp,*psp;
  signed char *orp;
  int *tsp;
  unsigned int i,j,k,l,m,n,p,nr,nl,ir,jr,n_ligs,n_occupied;
  int p_index,rx_index,subvol;
  int grid_size,uu,vv,vv_max;
  byte grid_shape;
  double area;
  double diag_x,diag_y,r_slope,width,u_width;
  double u_factor,u_factor_2,v_factor,v_val,binding_factor;

  no_printf("Creating new effector grid...\n");

  vectorize(wp->vert[0],wp->vert[1],&ab);
  vectorize(wp->vert[1],wp->vert[2],&bc);
  if (wp->wall_shape==RECT_POLY) {
    grid_shape=RECTANGULAR;
    area=wp->area;
    no_printf("Area = %g\n",area);
    grid_size=(int) ceil(sqrt(area/2));
    if (grid_size==0 || grid_size>46340) {
      fprintf(log_file,"MCell: too many tiles in effector grid: %d\n",grid_size);
      fflush(log_file);
      return(1);
    }
    n=2*grid_size*grid_size;
    binding_factor=n/area;
    printf("binding_factor = %g\n",binding_factor);
    r_slope=(wp->length_first/wp->length_last);
    diag_x=wp->length_first;
    diag_y=wp->length_last;
    u_factor=grid_size/wp->length_first;
    u_factor_2=wp->length_first/grid_size;
    v_factor=grid_size/wp->length_last;
    u_axis.x=1;
    u_axis.y=0;
    u_axis.z=0;
    v_axis.x=0;
    v_axis.y=1;
    v_axis.z=0;
    vectorize(wp->vert[0],wp->vert[1],&i_axis);
    normalize(&i_axis);
    vectorize(wp->vert[0],wp->vert[3],&j_axis);
    normalize(&j_axis);
  }
  else {
    grid_shape=TRIANGULAR;
    vectorize(wp->vert[0],wp->vert[2],&ac);
    area=wp->area;
    no_printf("Area = %g\n",area);
    grid_size=(int) ceil(sqrt(area));
    if (grid_size==0 || grid_size>65536) {
      fprintf(log_file,"MCell: too many tiles in effector grid: %d\n",grid_size);
      fflush(log_file);
      return(1);
    }
    n=grid_size*grid_size;
    binding_factor=n/area;
    printf("binding_factor = %g\n",binding_factor);
    p_b.x=wp->length_first;
    p_b.y=0;
    p_b.z=0;
    width=dot_prod(&ab,&ac)/wp->length_first;
    p_c.x=width;
    p_c.y=2*area/wp->length_first;
    p_c.z=0;
    r_slope=(p_c.x/p_c.y);
    diag_x=p_c.x;
    diag_y=p_c.y;
    p_d.x=wp->vert[0]->x+(width*ab.x/wp->length_first);
    p_d.y=wp->vert[0]->y+(width*ab.y/wp->length_first);
    p_d.z=wp->vert[0]->z+(width*ab.z/wp->length_first);
    vectorize(wp->vert[0],wp->vert[1],&i_axis);
    normalize(&i_axis);
    vectorize(&p_d,wp->vert[2],&j_axis);
    normalize(&j_axis);
    u_axis.x=p_c.y-p_b.y;
    u_axis.y=p_b.x-p_c.x;
    u_axis.z=0;
    normalize(&u_axis);
    u_width=dot_prod(&p_b,&u_axis);
    u_factor=grid_size/u_width;
    u_factor_2=wp->length_first/grid_size;
    v_factor=grid_size/p_c.y;
    vectorize(&p_b,&p_c,&v_axis);
    v_axis.x=v_axis.x/grid_size;
    v_axis.y=v_axis.y/grid_size;
    v_axis.z=0;

  }
  step_u.x=ab.x/grid_size;
  step_u.y=ab.y/grid_size;
  step_u.z=ab.z/grid_size;
  step_v.x=bc.x/grid_size;
  step_v.y=bc.y/grid_size;
  step_v.z=bc.z/grid_size;
  vectorize(&step_u,&step_v,&diagonal);

  no_printf("Initializing effector grid size %d ...\n",n);
  fflush(log_file);

  if ((tiles=(struct rx **)malloc(n*sizeof(struct rx *)))==NULL) {
    fprintf(log_file,"MCell: cannot store tiles for effector grid: %d\n",n);
    fflush(log_file);
    return(1);
  }
  if ((tsp=(int *)malloc(n*sizeof(int)))==NULL) {
    fprintf(log_file,"MCell: cannot store time stamp data for effector grid: %d\n",n);
    fflush(log_file);
    return(1);
  }
  if ((dsp=(unsigned short *)malloc(n*sizeof(unsigned short)))==NULL) {
    fprintf(log_file,"MCell: cannot store desired state data for effector grid: %d\n",n);
    fflush(log_file);
    return(1);
  }
  if ((orp=(signed char *)malloc(n*sizeof(signed char)))==NULL) {
    fprintf(log_file,"MCell: cannot store orientation data\n");
    fflush(log_file);
    return(1);
  }

  for (i=0;i<n;i++) {
    tiles[i]=NULL;
    tsp[i]=(-INT_MAX);
    dsp[i]=0;
    orp[i]=0;
  }

  no_printf("Done initializing effector grid size %d\n",n);
  fflush(log_file);

  if ((ep=(struct effector *)malloc(sizeof(struct effector)))==NULL) {
    return(1);
  }
  wp->effectors=ep;
  wp->effectors->grid_shape=grid_shape;
  wp->effectors->grid_size=grid_size;
  wp->effectors->u_axis.x=u_axis.x;
  wp->effectors->u_axis.y=u_axis.y;
  wp->effectors->v_axis.x=v_axis.x;
  wp->effectors->v_axis.y=v_axis.y;
  wp->effectors->i_axis.x=i_axis.x;
  wp->effectors->i_axis.y=i_axis.y;
  wp->effectors->i_axis.z=i_axis.z;
  wp->effectors->j_axis.x=j_axis.x;
  wp->effectors->j_axis.y=j_axis.y;
  wp->effectors->j_axis.z=j_axis.z;
  wp->effectors->r_slope=r_slope;
  wp->effectors->diag.x=diag_x;
  wp->effectors->diag.y=diag_y;
  wp->effectors->u_factor=u_factor;
  wp->effectors->u_factor_2=u_factor_2;
  wp->effectors->v_factor=v_factor;
  wp->effectors->r_u_factor=1/u_factor;
  wp->effectors->r_v_factor=1/v_factor;
  wp->effectors->binding_factor=binding_factor;
  wp->effectors->n_tiles=n;
  wp->effectors->n_occupied=0;
  
  wp->effectors->i=0; /* deprecated */
  wp->effectors->j=0; /* deprecated */

  wp->effectors->tiles=tiles;
  wp->effectors->time_stamp=tsp;
  wp->effectors->desired_state=dsp;
  wp->effectors->prev_state=NULL;
  wp->effectors->index=n_effector_grids++;
  wp->effectors->n_types=0;
  wp->effectors->orient=orp;
  wp->effectors->set=0;
  wp->effectors->wall=wp;
	
  no_printf("Done creating new effector grid.\n");

  return(0);
}


/**
 * Creates an array pointers to reference effector elements of the
 * wp linked list (walls) directly.
 * This effectively flattens the linked list.
 */
int init_effector_table(struct wall *wp)
{
  if (n_effector_grids>0) {
    if ((effector_table=(struct effector **)malloc
	 ((n_effector_grids)*sizeof(struct effector *)))==NULL) {
      return(1);
    }
    while (wp) {
      if (wp->effectors) {
	effector_table[wp->effectors->index]=wp->effectors;
       }
      wp=wp->next_wall;
    }
  }
  return(0);
}


/**
 * Creates an array pointers to reference elements of the
 * parent_rx linked list directly.
 * This effectively flattens the linked list.
 */
int init_rx_table(struct parent_rx *prxp)
{
  
  if (n_rx_types>0) {
    if ((rx_table=(struct parent_rx **)malloc
         (n_rx_types*sizeof(struct parent_rx *)))==NULL) {
      return(1);
    }
    while (prxp) {
      rx_table[prxp->rx_index]=prxp;
      prxp=prxp->next;
    }
  }
  return(0);
}

/**
 * Creates an array pointers to reference elements of the
 * release_event_queue linked list directly.
 * This effectively flattens the linked list.
 */
int init_release_event_table(struct release_event_queue *reqp)
{
  
  if (n_release_events>0) {
    if ((release_event_table=(struct release_event_queue **)malloc
         (n_release_events*sizeof(struct release_event_queue *)))==NULL) {
      return(1);
    }
    while (reqp) {
      release_event_table[reqp->index]=reqp;
      reqp=reqp->next;
    }
  }
  return(0);
}

/**
 * Frees up the memory used by the value of a symbol table entry
 * after we're done with it.
 * Currently deals with OBJ::POLY_OBJ::{BOX_POLY|ORDERED_POLY} types.
 */
void destroy_sym_value(struct sym_table *sym) {
  struct object *objp;
  struct polygon_object *pop;
  struct ordered_poly *opp;
  struct box_poly *bpp;
  struct vector3 *vect3;
  int i;

no_printf("Starting GARBAGE COLLECTING\n");
no_printf("sym: %s\n", sym->name);
  switch (sym->sym_type) {
	case OBJ:
      objp = (struct object *)sym->value;

      switch (objp->object_type) {
		case POLY_OBJ:
		  pop = (struct polygon_object *)objp->obj;
		  switch (pop->list_type) {
			case ORDERED_POLY:
              opp = (struct ordered_poly *)pop->polygon_data;
			  for (i=0;i<opp->n_verts;i++) {
				free(opp->vertex[i]);
			  }
			  free(opp->vertex);
			  if (opp->normal != NULL) {
				  for (i=0;i<opp->n_verts;i++) {
					free(opp->normal[i]);
				  }
				  free(opp->normal);
			  }
			  free(opp);
			break;

			case BOX_POLY:
			  bpp = (struct box_poly *)pop->polygon_data;
			  free(bpp->llf);
			  free(bpp->urb);
			  free(bpp);
			break;
		  }
		  free(pop);
		break;
	  }
	break;
  }

no_printf("Done GARBAGE COLLECTING\n");
fflush(stdout);
}




/* **************************************************************** */

#endif



/***************************************************************************
rel_expr_grab_obj:
  In: release expression
      place to allocate memory for temporary void_list
  Out: a linked list containing all the objects referred to in the
       release expression (including duplcates), or NULL if there are
       no such objects.
***************************************************************************/

/* Not the most efficient due to slow merging, but it works. */
struct void_list* rel_expr_grab_obj(struct release_evaluator *root,struct mem_helper *voidmem)
{
  struct void_list *vl = NULL;
  struct void_list *vr = NULL;
  
  if (root->left != NULL)
  {
    if (root->op&REXP_LEFT_REGION)
    {
      vl = mem_get(voidmem);
      if (vl==NULL) return NULL;
      vl->data = ((struct region*)(root->left))->parent;
      vl->next=NULL;
    }
    else vl = rel_expr_grab_obj(root->left,voidmem);
  }
  if (root->right != NULL)
  {
    if (root->op&REXP_RIGHT_REGION)
    {
      vr = mem_get(voidmem);
      if (vr==NULL) return NULL;
      vr->data = ((struct region*)(root->right))->parent;
      vr->next=NULL;
    }
    else vr = rel_expr_grab_obj(root->right,voidmem);
  }
  
  if (vl==NULL)
  {
    if (vr==NULL) return NULL;
    return vr;
  }
  else if (vr==NULL)
  {
    return vl;
  }
  else
  {
    struct void_list *vp;
    
    for (vp=vl;vp->next!=NULL;vp=vp->next) {}
    
    vp->next = vr;
    
    return vl;
  }
  return NULL;
}


/***************************************************************************
find_unique_rev_objects:
  In: release expression
      place to store the number of unique objects we find
  Out: an array of pointers to each object listed in the release
       expression (no duplicates), or NULL if out of memory.  The
       second argument is set to the length of the array.
***************************************************************************/

struct object** find_unique_rev_objects(struct release_evaluator *root,int *n)
{
  struct object **o_array;
  struct void_list *vp,*vq;
  struct mem_helper *voidmem;
  int i;
  
  voidmem = create_mem(1024,sizeof(struct void_list));
  
  vp = rel_expr_grab_obj(root,voidmem);
  if (vp==NULL) return NULL;
  
  vp = void_list_sort(vp);
  
  for (i=1,vq=vp ; vq!=NULL && vq->next!=NULL ; vq=vq->next , i++)
  {
    while (vq->data == vq->next->data)
    {
      vq->next = vq->next->next;
      if (vq->next==NULL) break;
    }
  }
  
  if (vq==NULL) i--;
  *n = i;
  
  o_array = (struct object**)malloc(i*sizeof(struct object*));
  if (o_array==NULL) return NULL;
  
  for (i=0,vq=vp ; vq!=NULL ; vq=vq->next,i++)
  {
    o_array[i] = (struct object*)vq->data;
  }
  
  delete_mem(voidmem);
  
  return o_array;
}


/***************************************************************************
eval_rel_region_expr:
  In: release expression for a 2D region release
      the number of distinct objects in the world listed in the expression
      array of pointers to each of those objects
      array of pointers to bit arrays specifying which walls of each
        object are included in this release
  Out: 0 on success, 1 on failure.  On success, the bit arrays are set
       so that they indicate which walls of each object are included in
       this release site.
***************************************************************************/

int eval_rel_region_expr(struct release_evaluator *expr,int n,struct object **objs,struct bit_array **result)
{
  int i;
  char bit_op;
  
  if (expr->left!=NULL)
  {
    if (expr->op&REXP_LEFT_REGION)
    {
      i = void_array_search((void**)objs,n,((struct region*)(expr->left))->parent);
      result[i] = duplicate_bit_array( ((struct region*)(expr->left))->membership );
      if (result[i]==NULL) return 1;
    }
    else
    {
      i = eval_rel_region_expr(expr->left,n,objs,result);
      if (i) return 1;
    }
    
    if (expr->right==NULL)
    {
      if (expr->op&REXP_NO_OP) return 0;
      else return 1;
    }
    
    if (expr->op&REXP_RIGHT_REGION)
    {
      i = void_array_search((void**)objs,n,((struct region*)(expr->right))->parent);
      if (result[i]==NULL)
      {
        result[i] = duplicate_bit_array( ((struct region*)(expr->right))->membership );
        if (result[i]==NULL) return 1;
      }
      else
      {
        if (expr->op&REXP_UNION) bit_op = '|';
        else if (expr->op&REXP_SUBTRACTION) bit_op = '-';
        else if (expr->op&REXP_INTERSECTION) bit_op = '&';
        else return 1;

        bit_operation(result[i],((struct region*)(expr->right))->membership,bit_op);
      }
    }
    else
    {
      struct bit_array *res2[n];
      for (i=0;i<n;i++) res2[i]=NULL;
      
      i = eval_rel_region_expr(expr->right,n,objs,result);
      if (i) return 1;
      
      for (i=0;i<n;i++)
      {
        if (res2[i]!=NULL)
        {
          if (result[i]==NULL) result[i] = res2[i];
          else
          {
            if (expr->op&REXP_UNION) bit_op = '|';
            else if (expr->op&REXP_SUBTRACTION) bit_op = '-';
            else if (expr->op&REXP_INTERSECTION) bit_op = '&';
            else return 1;
            
            bit_operation(result[i],res2[i],bit_op);
            free_bit_array(res2[i]);
          }
        }
      }
    }
  }
  else return 1;  /* Left should always have something! */
  
  return 0;
}


/***************************************************************************
init_rel_region_data_2d:
  In: release data for a release of 2D molecules onto a region
  Out: 0 on success, 1 on failure.  A summary of all potentially available
       space over all objects contained in the region expression is
       generated and stored in arrays (typically of length equal to the
       number of walls in the region expression).
***************************************************************************/

int init_rel_region_data_2d(struct release_region_data *rrd)
{
  int i,j,k;
  struct polygon_object *po;

  rrd->owners = find_unique_rev_objects(rrd->expression , &(rrd->n_objects));
  if (rrd->owners==NULL)
  {
    fprintf(world->err_file,"Error: cannot find any objects for region release\n");
    return 1;
  }
  
  rrd->in_release = (struct bit_array**)malloc(rrd->n_objects*sizeof(struct bit_array*));
  if (rrd->in_release==NULL)
  {
    fprintf(world->err_file,"Error: out of memory creating region lists for 2D region releases\n");
    return 1;
  }
  for (i=0;i<rrd->n_objects;i++) rrd->in_release[i]=NULL;
  
  i = eval_rel_region_expr(rrd->expression,rrd->n_objects,rrd->owners,rrd->in_release);
  if (i)
  {
    fprintf(world->err_file,"Error: could not evaluate region expression.\n");
    return 1;
  }
  for (i=0;i<rrd->n_objects;i++)
  {
    if (rrd->in_release[i]==NULL) 
    {
      if (rrd->owners[i]==NULL)
      {
	fprintf(world->err_file,"Object %d of %d in region expression was not found!\n",i+1,rrd->n_objects);
	return 1;
      }
      fprintf(world->err_file,"Error: could not generate region data on object %s\n(Out of memory?)\n",rrd->owners[i]->sym->name);
      return 1;
    }
  }
  
  rrd->walls_per_obj = (int*)malloc(rrd->n_objects*sizeof(int));
  if (rrd->walls_per_obj==NULL)
  {
    fprintf(world->err_file,"Error: out of memory creating wall counts for 2D region releases\n");
    return 1;
  }
  
  rrd->n_walls_included=0;
  for (i=0;i<rrd->n_objects;i++)
  {
    rrd->walls_per_obj[i] = count_bits(rrd->in_release[i]);
    rrd->n_walls_included += rrd->walls_per_obj[i];
  }
  
  rrd->cum_area_list = (double*)malloc(rrd->n_walls_included*sizeof(double));
  rrd->wall_index = (int*)malloc(rrd->n_walls_included*sizeof(int));
  rrd->obj_index = (int*)malloc(rrd->n_walls_included*sizeof(int));
  if (rrd->cum_area_list==NULL || rrd->wall_index==NULL || rrd->obj_index==NULL)
  {
    fprintf(world->err_file,"Error: out of memory creating area lists for 2D region releases\n");
    return 1;
  }
  
  j = 0;
  for (i=0;i<rrd->n_objects;i++)
  {
    k = rrd->owners[i]->object_type;
    if (k != POLY_OBJ && k != BOX_OBJ)
    {
      fprintf(world->err_file,"Error: found a region on something that isn't a box or polygon object?\n");
      return 1;
    }
    po = (struct polygon_object*)(rrd->owners[i]->contents);
    for (k=0;k<po->n_walls;k++)
    {
      if (get_bit(rrd->in_release[i],k))
      {
        rrd->cum_area_list[j] = rrd->owners[i]->wall_p[k]->area;
        rrd->obj_index[j] = i;
        rrd->wall_index[j] = k;
        j++;
      }
    }
  }
  
  for (i=1;i<rrd->n_walls_included;i++)
  {
    rrd->cum_area_list[i] += rrd->cum_area_list[i-1];
  }
  
  return 0;
}


/***************************************************************************
create_region_bbox:
  In: a region
  Out: pointer to a 2-element array contining the LLF and URB corners of
       a bounding box around the region, or NULL if out of memory.
***************************************************************************/

struct vector3* create_region_bbox(struct region *r)
{
  int i,j,k;
  struct vector3 *bbox;
  struct vector3 *v;
  
  bbox = (struct vector3*) malloc(2*sizeof(struct vector3));
  if (bbox==NULL) return NULL;
  
  j=0;
  for (i=0;i<r->membership->nbits;i++)
  {
    if (get_bit(r->membership,i))
    {
      if (!j)
      {
        bbox[0].x = bbox[1].x = r->parent->wall_p[i]->vert[0]->x;
        bbox[0].y = bbox[1].y = r->parent->wall_p[i]->vert[0]->y;
        bbox[0].z = bbox[1].z = r->parent->wall_p[i]->vert[0]->z;
      }
      for (k=0;k<3;k++)
      {
        v = r->parent->wall_p[i]->vert[k];
        if (bbox[0].x > v->x) bbox[0].x = v->x;
        else if (bbox[1].x < v->x) bbox[1].x = v->x;
        if (bbox[0].y > v->y) bbox[0].y = v->y;
        else if (bbox[1].y < v->y) bbox[1].y = v->y;
        if (bbox[0].z > v->z) bbox[0].z = v->z;
        else if (bbox[1].z < v->z) bbox[1].z = v->z;
      }
      j++;
    }
  }
  
  return bbox;
}


/***************************************************************************
eval_rel_region_bbox:
  In: release expression for a 3D region release
      place to store LLF corner of the bounding box for the release
      place to store URB corner
  Out: 0 on success, 1 on failure.  Bounding box is set based on release
       expression (based boolean intersection of bounding boxes for each
       region).  The function reports failure if any region is unclosed.
***************************************************************************/

int eval_rel_region_bbox(struct release_evaluator *expr,struct vector3 *llf,struct vector3 *urb)
{
  int i;
  struct region *r;
  
  if (expr->left!=NULL)
  {
    if (expr->op&REXP_LEFT_REGION)
    {
      r = (struct region*)(expr->left);
      if (r->manifold_flag==MANIFOLD_UNCHECKED)
      {
        if (is_manifold(r)) r->manifold_flag = IS_MANIFOLD;
        else
        {
          fprintf(world->log_file,"Error--cannot release a 3D molecule inside an unclosed region\n");
          return 1;
        }
      }
      
      if (r->bbox==NULL)
      {
        r->bbox = create_region_bbox(r);
        if (r->bbox==NULL) return 1;
      }
      
      llf->x = r->bbox[0].x;
      llf->y = r->bbox[0].y;
      llf->z = r->bbox[0].z;
      urb->x = r->bbox[1].x;
      urb->y = r->bbox[1].y;
      urb->z = r->bbox[1].z;
    }
    else
    {
      i = eval_rel_region_bbox(expr->left,llf,urb);
      if (i) return 1;
    }
    
    if (expr->right==NULL)
    {
      if (expr->op&REXP_NO_OP) return 0;
      else return 1;
    }
    
    if (expr->op&REXP_SUBTRACTION) return 0;
    else
    {
      struct vector3 llf2;
      struct vector3 urb2;
      
      if (expr->op&REXP_RIGHT_REGION)
      {
        r = (struct region*)(expr->right);
        if (r->manifold_flag==MANIFOLD_UNCHECKED)
        {
          if (is_manifold(r)) r->manifold_flag = IS_MANIFOLD;
          else
          {
            fprintf(world->log_file,"Error--cannot release a 3D molecule inside an unclosed region\n");
            return 1;
          }
        }
        
        if (r->bbox==NULL)
        {
          r->bbox = create_region_bbox(r);
          if (r->bbox==NULL) return 1;          
        }

        llf2.x = r->bbox[0].x;
        llf2.y = r->bbox[0].y;
        llf2.z = r->bbox[0].z;
        urb2.x = r->bbox[1].x;
        urb2.y = r->bbox[1].y;
        urb2.z = r->bbox[1].z;
      }
      else
      {
        i = eval_rel_region_bbox(expr->right,&llf2,&urb2);
        if (i) return 1;
      }
      
      if (expr->op&REXP_UNION)
      {
        if (llf->x > llf2.x) llf->x = llf2.x;
        if (llf->y > llf2.y) llf->y = llf2.y;
        if (llf->z > llf2.z) llf->z = llf2.z;
        if (urb->x < urb2.x) urb->x = urb2.x;
        if (urb->y < urb2.y) urb->y = urb2.y;
        if (urb->z < urb2.z) urb->z = urb2.z;
      }
      else if (expr->op&REXP_INTERSECTION)
      {
        if (llf->x < llf2.x) llf->x = llf2.x;
        if (llf->y < llf2.y) llf->y = llf2.y;
        if (llf->z < llf2.z) llf->z = llf2.z;
        if (urb->x > urb2.x) urb->x = urb2.x;
        if (urb->y > urb2.y) urb->y = urb2.y;
        if (urb->z > urb2.z) urb->z = urb2.z;
      }
      else return 1;
    }
  }
  else return 1;  /* Left should always have something! */
  
  return 0;
}


/***************************************************************************
init_rel_region_data_3d:
  In: release region data structure
  Out: 0 on success, 1 on failure, -1 if there is no volume contained
       in the release (user error).  The release must be for a 3D release
       of molecules.  eval_rel_region_bbox is called to perform the
       initialization.
***************************************************************************/

int init_rel_region_data_3d(struct release_region_data *rrd)
{
  int i;
  
  rrd->n_walls_included = 0;
  
  i = eval_rel_region_bbox(rrd->expression,&(rrd->llf),&(rrd->urb));

  if (i) return 1;

  if (rrd->llf.x >= rrd->urb.x ||
      rrd->llf.y >= rrd->urb.y ||
      rrd->llf.z >= rrd->urb.z)
  {
    return -1;  /* Special signal to print out "nothing in here" error msg */
  }
  
  return 0;
}


/***************************************************************************
output_regrel_eval_tree:
  In: file to print tree on
      prefix string to put before the current branch of the tree
      prefix character for left half of current branch
      prefix character for right half of current branch
      release expression
  Out: no return value.  The tree is printed to the file.
***************************************************************************/

void output_relreg_eval_tree(FILE *f,char *prefix,char cA,char cB,struct release_evaluator *expr)
{
  int l = strlen(prefix);
  char my_op;
  
  if (expr->op&REXP_NO_OP)
  {
    fprintf(f,"%s >%s\n",prefix,((struct region*)(expr->left))->sym->name);
  }
  else
  {
    char prefixA[l+3];
    char prefixB[l+3];
    strncpy(prefixA,prefix,l);
    strncpy(prefixB,prefix,l);
    prefixA[l] = cA;
    prefixB[l] = cB;
    prefixA[l+1] = prefixB[l+1] = ' ';
    prefixA[l+2] = prefixB[l+2] = 0;
    
    if (expr->op&REXP_LEFT_REGION)
    {
      fprintf(f,"%s >%s\n",prefix,((struct region*)(expr->left))->sym->name);
    }
    else
    {
      output_relreg_eval_tree(f,prefixA,' ','|',expr->left);
    }
    
    my_op = '?';
    if (expr->op & REXP_UNION) my_op = '+';
    else if (expr->op & REXP_INTERSECTION) my_op = '*';
    else if (expr->op & REXP_SUBTRACTION) my_op = '-';
    
    fprintf(f,"%s%c\n",prefix,my_op);
    
    if (expr->op&REXP_RIGHT_REGION)
    {
      fprintf(f,"%s >%s\n",prefix,((struct region*)(expr->left))->sym->name);
    }
    else
    {
      output_relreg_eval_tree(f,prefixA,'|',' ',expr->right);
    }
  }
}
  

/***************************************************************************
init_releases:
  In: nothing
  Out: 0 on success, 1 on failure.  All release sites are initialized.
       Right now, the only release sites that need to be initialized are
       releases on regions.
***************************************************************************/

int init_releases()
{
  struct release_event_queue *req;
  struct abstract_element *ae;
  struct schedule_helper *sh;
  int i,j;
  
  for (sh=world->releaser ; sh!=NULL ; sh=sh->next_scale)
  {
    for (i=-1;i<sh->buf_len;i++)
    {
      for ( ae = (i==-1)?sh->current:sh->circ_buf_head[i] ; ae!=NULL ; ae=ae->next )
      {
        req = (struct release_event_queue*)ae;
        if (req->release_site->release_shape == SHAPE_REGION)
        {
          if ((req->release_site->mol_type->flags & NOT_FREE) == 0)
          {
            j = init_rel_region_data_3d(req->release_site->region_data);
            if (j==-1)
            {
              fprintf(world->err_file,"Region release site is empty!  Ignoring!  Evaluation tree:\n");
              output_relreg_eval_tree(world->err_file," ",' ',' ',req->release_site->region_data->expression);
              req->release_site->release_number_method=CONSTNUM;
              req->release_site->release_number=0;
            }
            else if (j)
	    {
	      fprintf(world->err_file,"Error initializing region release for molecule %s\nEvaluation tree (3D release):\n",req->release_site->mol_type->sym->name);
	      output_relreg_eval_tree(world->err_file," ",' ',' ',req->release_site->region_data->expression);
	      return 1;
	    }
          }
          else
          {
            j = init_rel_region_data_2d(req->release_site->region_data);
            if (j)
	    { 
	      fprintf(world->err_file,"Error initializing region release for molecule %s\nEvaluation tree (2D release):\n",req->release_site->mol_type->sym->name);
	      output_relreg_eval_tree(world->err_file," ",' ',' ',req->release_site->region_data->expression);
	      return 1;
	    }
          }
        }
      }
    }
  }
  
  return 0;
}
