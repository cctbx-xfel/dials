from __future__ import division

class Test(object):

  def __init__(self):
    pass

  def run(self):
    self.tst_crystal()

  def tst_crystal(self):
    from dials.model.serialize.crystal import crystal_to_dict
    from dials.model.serialize.crystal import crystal_from_dict
    from dials.model.experiment import Crystal
    from scitbx import matrix
    real_space_a = matrix.col((35.2402102454, -7.60002142787, 22.080026774))
    real_space_b = matrix.col((22.659572494, 1.47163505925, -35.6586361881))
    real_space_c = matrix.col((5.29417246554, 38.9981792999, 4.97368666613))

    c1 = Crystal(
        real_space_a=real_space_a,
        real_space_b=real_space_b,
        real_space_c=real_space_c,
        space_group_symbol="P 1 2/m 1",
        mosaicity=0.1)

    d = crystal_to_dict(c1)
    c2 = crystal_from_dict(d)
    eps = 1e-7
    assert(abs(matrix.col(d['real_space_a']) - real_space_a) <= eps)
    assert(abs(matrix.col(d['real_space_b']) - real_space_b) <= eps)
    assert(abs(matrix.col(d['real_space_c']) - real_space_c) <= eps)
    assert(d['space_group_hall_symbol'] == "-P 2y")
    assert(d['mosaicity'] == 0.1)
    assert(c1 == c2)
    print 'OK'


if __name__ == '__main__':
  test = Test()
  test.run()
