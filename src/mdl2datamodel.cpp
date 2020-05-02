/*
 * mdl2datamodel.cpp
 *
 *  Created on: Mar 16, 2020
 *      Author: ahusar
 */
#include <string>
#include <iostream>
#include <fstream>
#include <assert.h>

#include "mdl2datamodel.h"
#include "datamodel_defines.h"
#include "json/json.h"

using namespace std;
using namespace Json;


// ------------- utilities -----------------------
static const char* get_sym_name(const sym_entry *s) {
  if (s == nullptr) {
    return "NULL SYMBOL";
  }
  else {
    return s->name;
  }
}

static bool has_flag(unsigned int flags, unsigned int one_flag) {
  assert(__builtin_popcount(one_flag) == 1); // check that we are testing just a single bit
  return (flags & one_flag) != 0;
}


// ------------- conversion methods -----------------------

void Mdl2DataModel::add_model_language(Value& mcell, const char* ver) {
  mcell[KEY_MODEL_LANGUAGE] = ver;
}

void Mdl2DataModel::api_add_version(Value& mcell, const int ver) {
  mcell[KEY_API_VERSION] = ver;
}

void Mdl2DataModel::cellblender_add_version(Value& mcell, const char* ver) {
  mcell[KEY_CELLBLENDER_VERSION] = ver;
}

void Mdl2DataModel::cellblender_add_sha1(Value& mcell, const char* ver) {
  mcell[KEY_CELLBLENDER_SHA1] = ver;
}

void Mdl2DataModel::blender_add_version(Value& mcell) {
  Value& blender_version = mcell[KEY_BLENDER_VERSION];
  if (BLENDER_VERSION == nullptr) return;
  blender_version.append(BLENDER_VERSION[0]);
  blender_version.append(BLENDER_VERSION[1]);
  blender_version.append(BLENDER_VERSION[2]);
}

bool Mdl2DataModel::convert_molecule(Value& molecule_list, species* spec) {

  const char* name = get_sym_name(spec->sym);
  if (strcmp(name, ALL_MOLECULES) == 0 ||
      strcmp(name, ALL_VOLUME_MOLECULES) == 0 ||
      strcmp(name, ALL_SURFACE_MOLECULES) == 0) {
    return true;
  }

  Value m;
  m[KEY_MOL_NAME] = name;
  m[KEY_DIFFUSION_CONSTANT] = spec->D;
  m[KEY_MOL_TYPE] = (has_flag(spec->flags, ON_GRID)) ? VALUE_MOL_TYPE_2D : VALUE_MOL_TYPE_3D;

  molecule_list.append(m);

  return true;
}

bool Mdl2DataModel::convert_molecule_list(Value& mcell) {
  Value& define_molecules = mcell[KEY_DEFINE_MOLECULES];
  Value& molecule_list = define_molecules[KEY_MOLECULE_LIST];

  for (int i = 0; i < s->n_species; i++) {
    species* spec = s->species_list[i];
    convert_molecule(molecule_list, spec);
  }

  json_add_version(define_molecules, JSON_DM_VERSION_1638);

  return true;
}

bool Mdl2DataModel::convert_reaction(Value& reaction_list, rxn *rx) {
  Value r;
  r[KEY_DATA_MODEL_VERSION] = JSON_DM_VERSION_1330; // hardcoded for now
  r[KEY_RXN_NAME] = rx->pathway_head->pathname->sym->name; // Is this the correct name? There are multiple options.

  // converting double to scientific notation string
  ostringstream streamObj;
  streamObj << rx->pathway_head->km;
  string fwd_rate = streamObj.str();
  r[KEY_RXN_FWD_RATE] = fwd_rate;

  reaction_list.append(r);

  return true;
}

bool Mdl2DataModel::convert_reaction_list(Value& mcell) {
  Value& define_reactions = mcell[KEY_DEFINE_REACTIONS];
  Value& reaction_list = define_reactions[KEY_REACTION_LIST];

  rxn** reaction_hash = s->reaction_hash;
  int count = s->rx_hashsize;

  for (int i = 0; i < count; ++i) {
    rxn *rxn_curr = reaction_hash[i];
    while (rxn_curr != nullptr) {
      convert_reaction(reaction_list, rxn_curr);
      rxn_curr = rxn_curr->next;
    }
  }

  json_add_version(define_reactions, JSON_DM_VERSION_1330);

  return true;
}

bool Mdl2DataModel::convert_object(Value& object_list, Value& vertex_list, Value& elem_conn_list, geom_object *obj) {
  if (obj == nullptr) {
    return true;
  }

  geom_object *obj_f = obj->first_child;
  double vert_list[obj_f->n_verts][3];

  // Appending vertex_list
  for(int i=0; i < obj_f->n_verts; i++) {
    vert_list[i][0] = obj_f->vertices[i]->x * VECTOR_SCALE_FACTOR;
    vert_list[i][1] = obj_f->vertices[i]->y * VECTOR_SCALE_FACTOR;
    vert_list[i][2] = obj_f->vertices[i]->z * VECTOR_SCALE_FACTOR;

    // inserting whole arrays into vertex_list
    Value tmp;
    for(auto j : vert_list[i]) { tmp.append(j); }
    vertex_list.append(tmp);
  }

  //Value o;
  //o[KEY_NAME] = *obj->last_name;
  //object_list.append(o); // this is throwing an instance of 'Json::LogicError'

  return true;
}

bool Mdl2DataModel::convert_geometrical_objects(Value& mcell) {
  Value& geometrical_objects = mcell[KEY_GEOMETRICAL_OBJECTS];
  Value& object_list = geometrical_objects[KEY_OBJECT_LIST];
  Value& vertex_list = object_list[KEY_VERTEX_LIST];
  Value& elem_conn_list = object_list[KEY_ELEMENT_CONNECTIONS];

  geom_object *root_inst = s->root_instance;
  for (geom_object* curr_obj = root_inst->first_child; curr_obj != nullptr; curr_obj = curr_obj->next) {
    convert_object(object_list, vertex_list, elem_conn_list, curr_obj);
  }

  return true;
}

bool Mdl2DataModel::convert(const volume* state_, const char* output_file_name) {

  s = state_;

  cout << "Convert called\n";

  Value root;
  Value& mcell = root[KEY_MCELL];

  // Calling datamodel conversion functions:
  cellblender_add_version(mcell, VALUE_CELLBLENDER_VERSION);
  convert_reaction_list(mcell);
  convert_geometrical_objects(mcell);
  add_model_language(mcell, VALUE_MCELL3);
  blender_add_version(mcell);
  convert_molecule_list(mcell);
  api_add_version(mcell, VALUE_API_VERSION_0);
  cellblender_add_sha1(mcell, VALUE_CELLBLENDER_SHA1);

  Json::StreamWriterBuilder wbuilder;
  wbuilder["indentation"] = " ";
  std::string document = Json::writeString(wbuilder, root);

  // write result into a file
  ofstream myfile(output_file_name);
  if (myfile.is_open())
  {
    myfile << document;
    myfile.close();
  }
  else {
    cout << "Unable to open file " << output_file_name << " for writing.\n";
  }

  return true;
}

// C-language style entry point
bool convert_to_datamodel(const volume* state, const char* output_file_name) {
  Mdl2DataModel converter;
  return converter.convert(state, output_file_name);
}
