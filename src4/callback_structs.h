/*
 * callback_structs.h
 *
 *  Created on: Oct 14, 2019
 *      Author: ahusar
 */

#ifndef SRC4_CALLBACK_STRUCTS_H_
#define SRC4_CALLBACK_STRUCTS_H_

#ifndef SWIG
#include "defines.h"
#endif

namespace MCell {

struct VEC4_ALIGNMENT WallHitInfo {
  vec3_t pos;
  molecule_id_t molecule_id;
  geometry_object_id_t geometry_object_id;
  wall_id_t wall_id;
  float_t time;
};

} /* namespace MCell */

#endif /* SRC4_CALLBACK_STRUCTS_H_ */
