
#ifndef DIALS_MODEL_SERIALIZE_SHOEBOX_FILE_IMPORTER_H
#define DIALS_MODEL_SERIALIZE_SHOEBOX_FILE_IMPORTER_H

#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/model/serialize/shoebox_file.h>
#include <dials/model/data/shoebox.h>

namespace dials { namespace model { namespace serialize {

  /**
   * A class to import shoeboxes
   */
  class ShoeboxFileImporter {
  public:

    typedef af::versa< double, af::c_grid<2> > gain_map_type;
    typedef af::versa< double, af::c_grid<2> > dark_map_type;
    typedef af::versa< bool, af::c_grid<2> > mask_map_type;
    
    typedef af::const_ref< double, af::c_grid<2> > gain_map_ref_type;
    typedef af::const_ref< double, af::c_grid<2> > dark_map_ref_type;
    typedef af::const_ref< bool, af::c_grid<2> > mask_map_ref_type;

    typedef af::versa< int, af::c_grid<3> > raw_shoebox_type;

    /**
     * Initialise the importer
     * @param filename The path to the shoebox file
     * @param gain The gain map
     * @param dark The dark map
     * @param mask The mask
     */
    ShoeboxFileImporter(const std::string &filename,
                        const af::const_ref<gain_map_ref_type> &gain,
                        const af::const_ref<dark_map_ref_type> &dark,
                        const af::const_ref<mask_map_ref_type> &mask)
        : reader_(filename),
          bbox_(reader_.bbox()),
          panel_(reader_.panels()) {

      // Ensure the number of maps match the number of panels
      DIALS_ASSERT(gain.size() > af::max(panel_.const_ref()));
      DIALS_ASSERT(gain.size() == dark.size());
      DIALS_ASSERT(gain.size() == mask.size());

      // Copy all the maps
      for (std::size_t i = 0; i < gain.size(); ++i) {

        // Ensure the sizes match
        DIALS_ASSERT(gain[i].accessor().all_gt(0));
        DIALS_ASSERT(gain[i].accessor().all_eq(dark[i].accessor()));
        DIALS_ASSERT(gain[i].accessor().all_eq(mask[i].accessor()));

        // Copy the data
        gain_map_type g(af::c_grid<2>(gain[i].accessor()));
        dark_map_type d(af::c_grid<2>(dark[i].accessor()));
        mask_map_type m(af::c_grid<2>(mask[i].accessor()));
        std::copy(gain[i].begin(), gain[i].end(), g.begin());
        std::copy(dark[i].begin(), dark[i].end(), d.begin());
        std::copy(mask[i].begin(), mask[i].end(), m.begin());

        // Append into arrays
        gain_maps_.push_back(g);
        dark_maps_.push_back(d);
        mask_maps_.push_back(m);
      }

      // Check bounding boxes do not go outside map ranges
      for (std::size_t i = 0; i < bbox_.size(); ++i) {
        gain_map_type::accessor_type acc = gain_maps_[panel_[i]].accessor();
        int6 b = bbox_[i];
        DIALS_ASSERT(b[0] >= 0 && b[2] >= 0 && b[4] >= 0);
        DIALS_ASSERT(b[1] <= acc[1] && b[3] <= acc[0]);
      }
    }

    /**
     * @returns The number of shoeboxes
     */
    std::size_t size() const {
      return reader_.size();
    }

    /**
     * @returns The bounding boxes
     */
    af::shared<int6> bboxes() const {
      return reader_.bbox();
    }

    /**
     * @returns The panel numbers
     */
    af::shared<std::size_t> panels() const {
      return reader_.panels();
    }

    /**
     * Read a shoebox
     * @param index The index of the shoebox
     * @returns The shoebox
     */
    Shoebox<> operator[](std::size_t index) {
      af::versa< int, af::c_grid<3> > data = reader_.read(index);

      // Ensure index is valid
      DIALS_ASSERT(index < size());

      // Get the panel and bounidng box
      std::size_t panel = panel_[index];
      int6 bbox = bbox_[index];

      // Initialise the shoebox
      Shoebox<> sbox(panel, bbox);

      // Allocate the memory
      sbox.allocate();

      // Grab the lookup tables
      gain_map_type gain = gain_maps_[panel];
      dark_map_type dark = dark_maps_[panel];
      mask_map_type mask = mask_maps_[panel];

      // Assign the values of the mask and pixels
      raw_shoebox_type raw = reader_.read(index);

      // Assign all the pixel values
      std::size_t j0 = bbox[2];
      std::size_t i0 = bbox[0];
      std::size_t zsize = bbox[5] - bbox[4];
      std::size_t ysize = bbox[3] - bbox[2];
      std::size_t xsize = bbox[1] - bbox[0];
      for (std::size_t k = 0; k < zsize; ++k) {
        for (std::size_t j = 0; j < ysize; ++j) {
          for (std::size_t i = 0; i < xsize; ++i) {
            std::size_t j1 = j0 + j;
            std::size_t i1 = i0 + i;
            sbox.data(k,j,i) = gain(j1,i1) * raw(k,j,i) - dark(j1,i1);
            if (mask(j1,i1)) {
              sbox.mask(k,j,i) = Valid;
            } else {
              sbox.mask(k,j,i) = 0;
            }
          }
        }
      }

      // return the shoebox
      return sbox;
    }

    /**
     * Get a range of shoeboxes.
     * @param i0 The start of the range.
     * @param i1 The end of the range
     * @returns The list of shoeboxes
     */
    af::shared< Shoebox<> > select(std::size_t i0, std::size_t i1) {
      af::shared< Shoebox<> > result;
      for (std::size_t i = i0; i < i1; ++i) {
        result.push_back((*this)[i]);
      }
      return result;
    }

    /**
     * Get a selection of shoeboxes
     * @param index The shoebox indices
     * @returns The list of shoeboxes
     */
    af::shared< Shoebox<> > select(const af::const_ref<std::size_t> &index) {
      af::shared< Shoebox<> > result;
      for (std::size_t i = 0; i < index.size(); ++i) {
        result.push_back((*this)[index[i]]);
      }
      return result;
    }

  private:

    ShoeboxFileReader reader_;
    af::shared<int6> bbox_;
    af::shared<std::size_t> panel_;
    af::shared<gain_map_type> gain_maps_;
    af::shared<dark_map_type> dark_maps_;
    af::shared<mask_map_type> mask_maps_;
  };

}}} // namespace dials::model::serialize

#endif // DIALS_MODEL_SERIALIZE_SHOEBOX_FILE_IMPORTER_H
