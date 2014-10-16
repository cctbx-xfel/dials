
from __future__ import division
from dials.algorithms.profile_model.gaussian_rs.transform import rebin_pixels

class Test(object):

  def __init__(self):
    pass

  def run(self):
    self.tst_identical()
    self.tst_known_offset()
    self.tst_larger_input()
    self.tst_larger_output()
    self.tst_known_orientation()
    self.tst_conservation_of_counts()

  def tst_identical(self):
    from scitbx.array_family import flex
    from random import uniform

    # Set the size of the grid
    height = 10
    width = 10

    # Create the grid data
    grid = flex.double([uniform(0, 100) for i in range(height * width)])
    grid.reshape(flex.grid(height, width))

    # Create the grid coordinates
    xy = []
    for j in range(height + 1):
      for i in range(width + 1):
        xy.append((i, j))
    gridxy = flex.vec2_double(xy)
    gridxy.reshape(flex.grid(height+1, width+1))

    # Get the output grid
    output = rebin_pixels(grid, gridxy, (height, width))

    # Check that each pixel is identical
    eps = 1e-7
    for j in range(height):
      for i in range(width):
        assert(abs(output[j,i] - grid[j, i]) <= eps)

    # Test passed
    print 'OK'

  def tst_known_offset(self):
    from scitbx.array_family import flex

    # Set the size of the grid
    height = 10
    width = 10

    # Create the grid data
    grid = flex.double([1 for i in range(height * width)])
    grid.reshape(flex.grid(height, width))

    # Create the grid coordinates
    xy = []
    for j in range(height + 1):
      for i in range(width + 1):
        xy.append((i + 0.5, j + 0.5))
    gridxy = flex.vec2_double(xy)
    gridxy.reshape(flex.grid(height+1, width+1))

    # Get the output grid
    output = rebin_pixels(grid, gridxy, (height, width))

    # Check that each each pixel along the left and bottom has a value
    # of 0.5 and that everything else is 1
    eps = 1e-7
    assert(abs(output[0,0] - 0.25) <= eps)
    for i in range(1, width):
      assert(abs(output[0,i] - 0.5) <= eps)
    for j in range(1, height):
      assert(abs(output[j,0] - 0.5) <= eps)
    for j in range(1, height):
      for i in range(1, width):
        assert(abs(output[j,i] - 1.0) <= eps)

    # Test passed
    print 'OK'

  def tst_larger_output(self):
    from scitbx.array_family import flex

    # Set the size of the grid
    height = 10
    width = 10

    # Create the grid data
    value = 13
    grid = flex.double([value for i in range(height * width)])
    grid.reshape(flex.grid(height, width))

    # Create the grid coordinates
    xy = []
    for j in range(height + 1):
      for i in range(width + 1):
        xy.append((i*2, j*2))
    gridxy = flex.vec2_double(xy)
    gridxy.reshape(flex.grid(height+1, width+1))

    # Get the output grid
    output = rebin_pixels(grid, gridxy, (height*2, width*2))

    # Check that each each pixel has a value of 0.25 the input
    eps = 1e-7
    for j in range(1, height):
      for i in range(1, width):
        assert(abs(output[j,i] - 0.25 * value) <= eps)

    # Test passed
    print 'OK'

  def tst_larger_input(self):
    from scitbx.array_family import flex

    # Set the size of the grid
    input_height = 20
    input_width = 20
    output_height = 10
    output_width = 10

    # Create the grid data
    value = 13
    grid = flex.double([value for i in range(input_height * input_width)])
    grid.reshape(flex.grid(input_height, input_width))

    # Create the grid coordinates
    xy = []
    for j in range(input_height + 1):
      for i in range(input_width + 1):
        xy.append((i/2.0, j/2.0))
    gridxy = flex.vec2_double(xy)
    gridxy.reshape(flex.grid(input_height+1, input_width+1))

    # Get the output grid
    output = rebin_pixels(grid, gridxy, (output_height, output_width))

    # Check that each each pixel has a value of 4 times the input
    eps = 1e-7
    for j in range(1, output_height):
      for i in range(1, output_width):
        assert(abs(output[j,i] - 4 * value) <= eps)

    # Test passed
    print 'OK'

  def tst_known_orientation(self):
    from scitbx.array_family import flex
    from scitbx import matrix
    from math import sin, cos, pi, sqrt

    # Set the size of the grid
    input_height = 4
    input_width = 4
    output_height = 4
    output_width = 4

    # Create the grid data
    value = 13
    grid = flex.double([value for i in range(input_height * input_width)])
    grid.reshape(flex.grid(input_height, input_width))

    # Create the grid coordinates
    xy = []
    R = matrix.sqr((cos(pi / 4), -sin(pi / 4), sin(pi / 4), cos(pi / 4)))
    for j in range(input_height + 1):
      for i in range(input_width + 1):
        ij = R * matrix.col((i * sqrt(8) / 4, j * sqrt(8) / 4))
        xy.append((ij[0] + 2, ij[1]))
    gridxy = flex.vec2_double(xy)
    gridxy.reshape(flex.grid(input_height+1, input_width+1))

    # Get the output grid
    output = rebin_pixels(grid, gridxy, (output_height, output_width))

    expected = [
      [0, 1, 1, 0],
      [1, 2, 2, 1],
      [1, 2, 2, 1],
      [0, 1, 1, 0]]

    # Check that each each pixel has a value of the input
    eps = 1e-7
    for j in range(1, output_height):
      for i in range(1, output_width):
        assert(abs(output[j,i] - expected[j][i] * value) <= eps)

    # Test passed
    print 'OK'

  def tst_conservation_of_counts(self):
    from scitbx.array_family import flex
    from scitbx import matrix
    from math import sin, cos, pi
    from random import uniform

    # Set the size of the grid
    input_height = 10
    input_width = 10
    output_height = 50
    output_width = 50

    # Create the grid data
    grid = flex.double([uniform(0, 100)
        for i in range(input_height * input_width)])
    grid.reshape(flex.grid(input_height, input_width))

    # Create the grid coordinates
    xy = []
    angle = uniform(0, pi)
    offset = (uniform(20, 30), uniform(20, 30))
    R = matrix.sqr((cos(angle), -sin(angle), sin(angle), cos(angle)))
    for j in range(input_height + 1):
      for i in range(input_width + 1):
        ij = R * matrix.col((i, j))
        xy.append((ij[0] + offset[0], ij[1] + offset[0]))
    gridxy = flex.vec2_double(xy)
    gridxy.reshape(flex.grid(input_height+1, input_width+1))

    # Get the output grid
    output = rebin_pixels(grid, gridxy, (output_height, output_width))

    # Check that the sum of the counts is conserved
    eps = 1e-7
    assert(abs(flex.sum(output) - flex.sum(grid)) <= eps)

    # Test passed
    print 'OK'

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
