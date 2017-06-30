from __future__ import division

from libtbx import easy_pickle
from dials.array_family import flex
from scitbx import matrix
import math, dxtbx

class PhiMask(object):
  def __init__(self, image_path):
    self.image = dxtbx.load(image_path)
    self.x = matrix.col((1, 0, 0))
    self.y = matrix.col((0, 1, 0))
    self.beam = -matrix.col(self.image.get_beam().get_s0()).normalize()
    assert self.beam == matrix.col((0, 0, 1))
    self.get_pixel_phi_array()

  def get_pixel_phi_array(self):
    self.all_phi = []
    for panel in self.image.get_detector():
      panel_fast, panel_slow = panel.get_image_size()
      size = panel_fast*panel_slow
      pix_range = flex.int(xrange(size))
      x_vec = flex.vec3_double(size, self.x)
      y_vec = flex.vec3_double(size, self.y)
      fast_coords = pix_range % panel_fast
      slow_coords = pix_range / panel_fast
      lab_coords = panel.get_lab_coord(panel.pixel_to_millimeter(
                     flex.vec2_double(fast_coords.as_double(), slow_coords.as_double())))
      x_projections = lab_coords.dot(x_vec)*x_vec
      y_projections = lab_coords.dot(y_vec)*y_vec
      sel = lab_coords.dot(y_vec) < 0
      xy_projections = x_projections + y_projections
      phi = xy_projections.angle(x_vec)
      phi.set_selected(sel, (2*math.pi)-phi.select(sel))
      phi.reshape(flex.grid((panel_fast, panel_slow)))
      self.all_phi.append(phi)

  def get_phi_mask(self, phi_start, phi_end):
    masks = []
    for phi_array in self.all_phi:
      keep = (phi_array >= phi_start) & (phi_array < phi_end)
      masks.append(keep)
    return masks

  def write_masks_by_phi(self, n_wedges, name="mask_by_phi"):
    phi_sequence_radians = [2*math.pi*i/n_wedges for i in xrange(n_wedges)]
    for i in xrange(len(phi_sequence_radians)):
      phi_start = phi_sequence_radians[i]
      phi_end = phi_start + phi_sequence_radians[1]
      masks = self.get_phi_mask(phi_start, phi_end)
      print "Saving mask between phi %.1f and %.1f"%(phi_start*180/math.pi, phi_end*180/math.pi)
      easy_pickle.dump("%s_%d.pickle" % (name, i), tuple(masks))

def run(args):
  file_path, n_wedges = args

  masker = PhiMask(file_path)
  masker.write_masks_by_phi(int(n_wedges))

if __name__ == "__main__":
  import sys
  run(sys.argv[1:])
