from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import unittest

import antenna

class TestEnhancedAntenna(unittest.TestCase):
  def test_standard_gain_ver_hor(self):
      #dirs = [{'hor':3.5,'ver':-30},{'hor':47.3,'ver':45},{'hor':342,'ver':70}]
      #gains = antenna.GetStandardAntennaGainsHorAndVer(dirs,30,10,peak_ant_gain = 5)
      #self.assertEqual(np.max(np.abs(gains - 5 * np.ones(4))), 0)
      dirs = {'hor':20,'ver':0}
      [hor_gain,ver_gain] = antenna.GetStandardAntennaGainsHorAndVer(dirs,ant_azimuth = 20,ant_mech_downtilt = 0,ant_hor_beamwidth = 120,ant_ver_beamwidth = 60,peak_ant_gain = 10)
      np.testing.assert_array_equal([hor_gain,ver_gain],[10,10])

      dirs = {'hor':60,'ver':30}
      [hor_gain,ver_gain] = antenna.GetStandardAntennaGainsHorAndVer(dirs,ant_azimuth = 0,ant_mech_downtilt = 0,ant_hor_beamwidth = 120,ant_ver_beamwidth = 60,peak_ant_gain = 10)
      np.testing.assert_array_equal([hor_gain,ver_gain],[10-3,10-3])

      dirs = {'hor':-15,'ver':9}
      [hor_gain,ver_gain] = antenna.GetStandardAntennaGainsHorAndVer(dirs,ant_azimuth = 45,ant_mech_downtilt = 2,ant_hor_beamwidth = 120,ant_ver_beamwidth = 20,peak_ant_gain = 10)
      np.testing.assert_array_equal([hor_gain,ver_gain],[10-3,10-3])

  def test_GetAntennaGainsFromGivenPattern(self):
      hor_pattern = {}
      hor_pattern['angle'] = list(range(-180,180,1))
      hor_pattern['gain'] = list(range(0,360,1))
      
      ver_pattern = {}
      ver_pattern['angle'] = list(range(-90,90,1))
      ver_pattern['gain'] = list(range(0,180,1))
      dirs = {'hor':20.5,'ver':10.5}
      [G_H_theta_R, G_V_phi_R, G_V_phi_Rsup] = antenna.GetAntennaGainsFromGivenPattern(dirs,hor_pattern,ver_pattern,ant_azimuth = 0,ant_mech_downtilt = 0)  
      np.testing.assert_array_equal([G_H_theta_R, G_V_phi_R, G_V_phi_Rsup],[200.5,100.5,100.5])
      
  def test_GetTwoDimensionalAntennaGain(self):
      hor_pattern = {}
      hor_pattern['angle'] = list(range(-180,180,1))
      hor_pattern['gain'] = list(range(0,360,1))
      
      ver_pattern = {}
      ver_pattern['angle'] = list(range(-90,90,1))
      ver_pattern['gain'] = list(range(0,180,1))
      dirs = {'hor':20.5,'ver':10.5}
      [G_H_theta_R, G_V_phi_R, G_V_phi_Rsup] = antenna.GetAntennaGainsFromGivenPattern(dirs,hor_pattern,ver_pattern,ant_azimuth = 00,ant_mech_downtilt = 0) 

      gain_two_dimensional = antenna.GetTwoDimensionalAntennaGain(dirs,G_H_theta_R,G_V_phi_R,G_V_phi_Rsup,
                                                                  hor_pattern,ver_pattern,ant_azimuth = None,
                                                                  ant_mech_downtilt=None,ant_elec_downtilt=0,peak_ant_gain=0)
       

if __name__ == '__main__':
  unittest.main()