/*
 * integration_ext.cc
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_xds_coordinate_system();
  void export_from_beam_vector_to_xds();
  void export_from_xds_to_beam_vector();
  void export_from_xds_e3_to_phi();
  void export_bbox_calculator();
  void export_find_overlapping_reflections();
  void export_reflection_profile_helpers();
  void export_shoebox_masker();
  void export_summation();
  void export_reciprocal_space_transform();

  BOOST_PYTHON_MODULE(dials_algorithms_integration_ext)
  {
    export_xds_coordinate_system();
    export_from_beam_vector_to_xds();
    export_from_xds_to_beam_vector();
    export_from_xds_e3_to_phi();
    export_bbox_calculator();
    export_find_overlapping_reflections();
    export_reflection_profile_helpers();
    export_shoebox_masker();
    export_summation();
    export_reciprocal_space_transform();
  }

}}} // namespace = dials::algorithms::boost_python
