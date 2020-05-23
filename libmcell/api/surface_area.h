/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
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

#ifndef API_SURFACE_REGION_H
#define API_SURFACE_REGION_H

#include "generated/gen_surface_area.h"
#include "api/common.h"
#include "defines.h"

namespace MCell {
namespace API {

class SurfaceArea: public GenSurfaceArea {
public:
  SURFACE_AREA_CTOR()

  void postprocess_in_ctor() {
    region_id = REGION_ID_INVALID;
  }

  std::shared_ptr<Region> as_region() override {
    assert("TODO" && false);
  }

  // simulation engine mapping
  region_id_t region_id;
};

} // namespace API
} // namespace MCell

#endif // API_SURFACE_REGION_H
