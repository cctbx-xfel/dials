# LIBTBX_SET_DISPATCHER_NAME dev.dials.semisynthetic_variance_analysis

from __future__ import division

def weighted_mean_variance(values, variances):
  import math
  weights = [1.0 / v for v in variances]
  sum_weights = sum(weights)
  weighted_mean = sum([v * w for v, w in zip(values, weights)]) / sum_weights
  weighted_variance = len(values) / sum_weights
  return weighted_mean, weighted_variance

def npp(values, input_mean_variance):
  import math
  from scitbx.math import distributions
  from scitbx.array_family import flex
  distribution = distributions.normal_distribution()
  values = flex.sorted(values)
  mean, variance = input_mean_variance
  scaled = (values - mean) / math.sqrt(variance)
  expected = distribution.quantiles(values.size())

  return expected, scaled

def semisynthetic_variance_analysis(semisynthetic_integrated_data_files):
  import cPickle as pickle
  from logging import info
  from dials.array_family import flex
  from dials.util.add_hash import add_hash, dehash

  integrated_data_sets = [pickle.load(open(data_file, 'rb')) for
                          data_file in semisynthetic_integrated_data_files]

  # first prepare the data files i.e. remove partials, keep only integrated
  # reflections, add the hash column, add weight column

  hash_set = None

  hashed_data_sets = []

  for integrated_data in integrated_data_sets:
    sel = integrated_data.get_flags(integrated_data.flags.integrated)
    integrated_data = integrated_data.select(sel)
    sel = integrated_data['partiality'] > 0.99
    integrated_data = integrated_data.select(sel)
    integrated_data = add_hash(integrated_data)
    hashed_data_sets.append(integrated_data)
    if hash_set is None:
      hash_set = set(integrated_data['hash'])
    else:
      hash_set = hash_set.intersection(set(integrated_data['hash']))

  # now analyse those reflections found to be in all data sets (here looking
  # at the profile fitted intensity and variance thereof)

  for h in hash_set:
    values_profile = flex.double()
    variances_profile = flex.double()
    for i in hashed_data_sets:
      sel = i['hash'] == h
      isel = sel.iselection()
      assert(len(isel) == 1)
      values_profile.append(i[isel[0]]['intensity.prf.value'])
      variances_profile.append(i[isel[0]]['intensity.prf.variance'])
    weighted_mean, weighted_variance = weighted_mean_variance(values_profile,
                                                              variances_profile)
    expected, scaled = npp(values_profile, (weighted_mean, weighted_variance))
    for e, s in zip(expected, scaled):
      print e, s

if __name__ == '__main__':
  import sys
  semisynthetic_variance_analysis(sys.argv[1:])