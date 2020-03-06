/******************************************************************************
 *
 * Copyright (C) 2019 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 *
******************************************************************************/

#ifndef SRC4_DYN_VERTEX_STRUCTS_H_
#define SRC4_DYN_VERTEX_STRUCTS_H_

#include "defines.h"

namespace MCell {


class Partition;

struct VEC4_ALIGNMENT VertexMoveInfo {
  VertexMoveInfo(const vertex_index_t vertex_index_, const vec3_t& translation_vec_)
    : translation_vec(translation_vec_), vertex_index(vertex_index_) {
  }
  // by how much
  vec3_t translation_vec;
  // and which index to move
  vertex_index_t vertex_index;
};

typedef std::vector<VertexMoveInfo, vec4_allocator<VertexMoveInfo>> VertexMoveInfoVector;
typedef std::map<wall_index_t, VertexMoveInfoVector> WallsWithTheirMovesMap;

struct VolumeMoleculeMoveInfo {
  VolumeMoleculeMoveInfo(const molecule_id_t molecule_id_, const wall_index_t wall_index_, const bool place_above_)
    : molecule_id(molecule_id_), wall_index(wall_index_), place_above(place_above_) {
  }
  // molecule to move
  molecule_id_t molecule_id;
  // which wall moved this molecule first
  wall_index_t wall_index;
  // above or below
  bool place_above;
};
typedef std::vector<VolumeMoleculeMoveInfo> VolumeMoleculeMoveInfoVector;

struct VEC4_ALIGNMENT SurfaceMoleculeMoveInfo {
  SurfaceMoleculeMoveInfo(const molecule_id_t molecule_id_, const wall_index_t wall_index_, const vec3_t pos3d_)
    : molecule_id(molecule_id_), wall_index(wall_index_), pos3d(pos3d_) {
  }
  // above or below
  vec3_t pos3d;
  // molecule to move
  molecule_id_t molecule_id;
  // which wall moved this molecule first
  wall_index_t wall_index;
};
typedef std::vector<SurfaceMoleculeMoveInfo, vec4_allocator<SurfaceMoleculeMoveInfo>> SurfaceMoleculeMoveInfoVector;

} // namespace MCell

#endif /* SRC4_DYN_VERTEX_STRUCTS_H_ */
