/*
 * mdl2datamodel.h
 *
 *  Created on: Mar 16, 2020
 *      Author: ahusar
 */

#ifndef SRC_MDL2DATAMODEL_H_
#define SRC_MDL2DATAMODEL_H_

#define VECTOR_SCALE_FACTOR 0.01

#include "mcell_structs.h"
#include "json/json.h"

class Mdl2DataModel {
public:
  bool convert(const volume* state_, const char* output_file_name);

private:
  bool convert_geometrical_objects(Json::Value& mcell);
  bool convert_object(Json::Value& object_list, Json::Value& vertex_list, Json::Value& elem_conn_list, geom_object *obj);
  bool convert_reaction_list(Json::Value& mcell);
  bool convert_reaction(Json::Value& reaction_list, rxn *rxn);
  bool convert_molecule(Json::Value& molecule_list, species* spec);
  bool convert_molecule_list(Json::Value& mcell);
  void cellblender_add_version(Json::Value& mcell, const char* ver);
  void blender_add_version(Json::Value& mcell);
  void add_model_language(Json::Value& mcell, const char* ver);
  void cellblender_add_sha1(Json::Value& mcell, const char* ver);
  void api_add_version(Json::Value& mcell, const int ver);

  // state initialized when convert
  const volume* s;
};

bool convert_to_datamodel(const volume* state, const char* output_file_name);

#endif /* SRC_MDL2DATAMODEL_H_ */
