#    Copyright 2017 SAS Project Authors. All Rights Reserved.
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.

"""Antenna gain routines for SAS.

Typical usage:
  # Get the gains for a given antenna horizontal pattern.
  # To be used for ESC and CBSD with actual antenna patterns.
  gains = GetAntennaPatternGains(hor_dirs, ant_azimuth, pattern)

  # Get CBSD antenna gain (from antenna defined by beamwidth and gain only).
  gains = GetStandardAntennaGains(hor_dirs, ant_azimuth, beamwidth, ant_gain)

  # Get Radar normalized antenna gain for DPA protection.
  gains = GetRadarNormalizedAntennaGains(hor_dirs, radar_azimuth)

  # Get FSS earth station antenna gains
  gains = GetFssAntennaGains(hor_dirs, ver_dirs,
                             fss_azimuth, fss_elevation, fss_ant_gain)
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np


def GetStandardAntennaGainsHorAndVer(dirs,ant_azimuth=None,ant_mech_downtilt=None,ant_elec_downtilt=0,ant_hor_beamwidth=None,ant_ver_beamwidth=None,
                                    ant_FBR=None,ant_USLS=None,peak_ant_gain=0):
  """Computes the antenna gain pattern using both horizontal and vertical properties.
  

  Directions and azimuth are defined compared to the north in clockwise
  direction and shall be within [0..360] degrees.

  Inputs:
    hor_dirs:       Ray directions in horizontal plane (degrees).
                    Either a scalar or an iterable.
    ver_dirs:       Ray directions in vertical plane (degrees).
                    Either a scalar or an iterable.
    ant_azimuth:     Antenna azimuth (degrees).
    ant_elevation:   Antenna elevation (degrees).
    ant_hor_beamwidth:  Antenna 3dB cutoff beamwidth (degrees) in horizontal plane
                    If None, then antenna is isotropic (default).
    ant_ver_beamwidth:  Antenna 3dB cutoff beamwidth (degrees) in vertical plane
                    If None, then antenna is isotropic (default).
    ant_mech_downtilt: Antenna mechanical downtilt(value withtin range 0 to +-15 degrees)
    ant_elec_downtilt: Antenna electrical downtilt
    peak_ant_gain:       Antenna gain (dBi)at boresight

  Returns:
    The CBSD antenna gains (in dBi).
    Either a scalar if hor_dirs and ver_dirs are scalar or an ndarray otherwise.
  """
  is_list = isinstance(dirs,list)
  hor_dirs = dirs['hor']
  ver_dirs = dirs['ver']
  hor_dirs = np.atleast_1d(hor_dirs)
  ver_dirs = np.atleast_1d(ver_dirs)
  if(ant_FBR is None):
    ant_FBR = 20
  if(ant_USLS is None):
    ant_USLS = 15

  if (ant_hor_beamwidth is None or ant_azimuth is None or
      ant_hor_beamwidth == 0 or ant_hor_beamwidth == 360):
    hor_gains = peak_ant_gain * np.ones(hor_dirs.shape)
  else:
    theta_R = hor_dirs - ant_azimuth
    theta_R[theta_R > 180] -= 360
    theta_R[theta_R < -180] += 360
    hor_gains = -min([12 * (theta_R / float(ant_hor_beamwidth))**2,ant_FBR])
    #gains[gains < -20] = -20.
    hor_gains += peak_ant_gain
  
  if (ant_ver_beamwidth is None or ant_mech_downtilt is None or
      ant_ver_beamwidth == 0 or ant_ver_beamwidth == 360):
    ver_gains = peak_ant_gain * np.ones(ver_dirs.shape)
  else:
    phi_R = ver_dirs + ant_mech_downtilt*np.cos(theta_R*np.pi/180) + ant_elec_downtilt
    
    ver_gains = -min([12 * (phi_R / float(ant_ver_beamwidth))**2,ant_USLS])
    #gains[gains < -20] = -20.
    ver_gains += peak_ant_gain

  if not is_list: return hor_gains[0],ver_gains[0]
  return hor_gains,ver_gains

def GetAntennaGainsFromGivenPattern(dirs,hor_pattern,ver_pattern, ant_azimuth = None,ant_mech_downtilt = None,ant_elec_downtilt = 0):


  """ REL2-R3-SGN-52105: Method B1 based Antenna Gain Calculation, step a


  Computes the gain at a given direction from a given antenna pattern(horizontal and vertical).
   
  Output gains will be calculated in horizontal and vertical dimension.

  Directions and azimuth are defined compared to the north in clockwise
  direction and shall be within [0..360] degrees.

  Inputs:
    hor_pattern: contains horizontal plane angles and associated gains 
    ver_pattern: contains horizontal plane angles and associated gains 
    hor_dir:  azimuth angle of the cbsd towards the receiver, relative to true north
    ver_dir:  elevation angle of the cbsd towrads the receiver
    ant_azimuth:     Antenna azimuth (degrees).
    ant_elevation:   Antenna elevation (degrees).

  Outputs:
    cbsd gain in horizontal and vertical plains at the given directions
  """
  
  is_list = isinstance(dirs,list)
  hor_dir = dirs['hor']
  ver_dir = dirs['ver']
  
  theta_R = hor_dir - ant_azimuth #azimuth angle of the line between the cbsd  main beam and receiver location, relative to cbsd antenna boresight
  if theta_R > 180:
    theta_R -= 360 
  elif theta_R < -180:
    theta_R += 360
  
  #elevation angle of the line between the cbsd  main beam and receiver location, relative to cbsd antenna boresight
  phi_R = ver_dir + ant_mech_downtilt*np.cos(theta_R*180/np.pi) + ant_elec_downtilt 
  
  theta_list = hor_pattern['angle']
  G_H_list = hor_pattern['gain']
  theta_R_idx = [i for i,j in enumerate(theta_list) if j == theta_R]

  phi_list = list(ver_pattern['angle'])
  phi_Rsup_list = list([180 -i for i in ver_pattern['angle']])
  
  G_V_list = list(ver_pattern['gain'])
  phi_R_idx = [i for i,j in enumerate(phi_list) if j == phi_R]  

  phi_R_supplementary_angle = 180 - phi_R
  phi_Rs_idx = [i for i,j in enumerate(phi_Rsup_list) if j == phi_R_supplementary_angle]  

  if any(theta_R_idx):
    G_H_theta_R = G_H_list[theta_R_idx]
  else:
    theta_diff = [theta_R - i for i in theta_list]
    theta_diff_pos = [i for i in theta_diff if i>0]
    theta_m = theta_list[theta_diff.index(min(theta_diff_pos))]

    theta_diff_neg = [i for i in theta_diff if i<0]
    
    theta_m_1 = theta_list[theta_diff.index(max(theta_diff_neg))]

    theta_m_idx = [i for i,j in enumerate(theta_list) if j == theta_m]
    G_H_theta_m = G_H_list[theta_m_idx[0]]

    theta_m_1_idx = [i for i,j in enumerate(theta_list) if j == theta_m_1]
    G_H_theta_m_1 = G_H_list[theta_m_1_idx[0]]

    G_H_theta_R_interp = ((theta_m_1-theta_R)*G_H_theta_m + (theta_R - theta_m)*G_H_theta_m_1)/(theta_m_1-theta_m)
    G_H_theta_R = G_H_theta_R_interp

  if any(phi_R_idx):
    G_V_phi_R = G_V_list[phi_R_idx]
  else:
    phi_diff = [phi_R - i for i in phi_list]
    phi_diff_pos = [i for i in phi_diff if i>0]
    phi_n = phi_list[phi_diff.index(min(phi_diff_pos))]

    phi_diff_neg = [i for i in phi_diff if i<0]
    
    phi_n_1 = phi_list[phi_diff.index(max(phi_diff_neg))]

    phi_n_idx = [i for i,j in enumerate(phi_list) if j == phi_n]
    G_V_phi_n = G_V_list[phi_n_idx[0]]

    phi_n_1_idx = [i for i,j in enumerate(phi_list) if j == phi_n_1]
    G_V_phi_n_1 = G_V_list[phi_n_1_idx[0]]

    G_V_phi_R_interp = ((phi_n_1-phi_R)*G_V_phi_n + (phi_R - phi_n)*G_V_phi_n_1)/(phi_n_1 - phi_n)
    G_V_phi_R = G_V_phi_R_interp

  if any(phi_Rs_idx):
    G_V_phi_Rsup = G_V_list[phi_Rs_idx]
  else:
    phi_Rsup_diff = [phi_R_supplementary_angle - i for i in phi_Rsup_list]
    phi_Rs_diff_pos = [i for i in phi_Rsup_diff if i>0]
    phi_k = phi_Rsup_list[phi_Rsup_diff.index(min(phi_Rs_diff_pos))]

    phi_Rs_diff_neg = [i for i in phi_Rsup_diff if i<0]
    phi_k_1 = phi_Rsup_list[phi_Rsup_diff.index(max(phi_Rs_diff_neg))]

    phi_k_idx = [i for i,j in enumerate(phi_Rsup_list) if j == phi_k]
    G_V_phi_k = G_V_list[phi_k_idx[0]]

    phi_k_1_idx = [i for i,j in enumerate(phi_Rsup_list) if j == phi_k_1]
    G_V_phi_k_1 = G_V_list[phi_k_1_idx[0]]

    G_V_phi_Rsup_interp = ((phi_k_1-phi_R_supplementary_angle)*G_V_phi_k + (phi_R_supplementary_angle - phi_k)*G_V_phi_k_1)/(phi_k_1 - phi_k)
    G_V_phi_Rsup = G_V_phi_Rsup_interp
  
    
  return G_H_theta_R, G_V_phi_R, G_V_phi_Rsup


def GetTwoDimensionalAntennaGain(hor_dir,ver_dir,hor_gain,ver_gain,ver_gain_sup_angle,
                                 hor_pattern,ver_pattern,ant_azimuth = None,ant_mech_downtilt=None,ant_elec_downtilt=0,peak_ant_gain=0):

                           
  """REL2-R3-SGN-52105: Method B1 based Antenna Gain Calculation, step b

  Computes the two dimensional antenna gain at a given direction, from horizontal and vertical gain.

  Directions and azimuth are defined compared to the north in clockwise
  direction and shall be within [0..360] degrees.

  Inputs:
    hor_dirs:       Ray directions in horizontal plane (degrees).
                    Either a scalar or an iterable.
    ver_dirs:       Ray directions in vertical plane (degrees).
                    Either a scalar or an iterable.
    ant_azimuth:     Antenna azimuth (degrees).
    ant_elevation:   Antenna elevation (degrees).
    ant_hor_beamwidth:  Antenna 3dB cutoff beamwidth (degrees) in horizontal plane
                    If None, then antenna is isotropic (default).
    ant_ver_beamwidth:  Antenna 3dB cutoff beamwidth (degrees) in vertical plane
                    If None, then antenna is isotropic (default).
    ant_mech_downtilt: Antenna mechanical downtilt(value withtin range 0 to +-15 degrees)
    ant_elec_downtilt: Antenna electrical downtilt
    peak_ant_gain:       Antenna gain (dBi)at boresight
  """
  G_H = hor_pattern["gain"]
  G_H_theta_R = hor_gain
  G_V = ver_pattern["gain"]
  G_V_phi_R = ver_gain
  G_V_phi_Rsup = ver_gain_sup_angle
  G_0 = peak_ant_gain
  G_cbsd_abs = G_H_theta_R + ( (1-abs(hor_dir)/180)*(G_V_phi_R - G_H[0]) + (abs(hor_gain)/180)*(G_V_phi_Rsup - G_H[179]))
  G_cbsd = G_cbsd_abs + G_0
  gain_two_dimensional = G_cbsd
  return gain_two_dimensional

def GetAntennaPatternGains(hor_dirs, ant_azimuth,
                           hor_pattern,
                           ant_gain=0):
  """Computes the gain for a given antenna pattern.

  Directions and azimuth are defined compared to the north in clockwise
  direction and shall be within [0..360] degrees.

  Inputs:
    hor_dirs:       Ray directions in horizontal plane (degrees).
                    Either a scalar or an iterable.
    ant_azimuth:    Antenna azimuth (degrees).
    hor_pattern:    The antenna horizontal pattern defined as a ndarray of
                    360 values with the following conventions:
                     - clockwise increments of 1 degree.
                     - 0 degree corresponds to the antenna boresight.
                     - values are 'gain' (and not attenuation)
    ant_gain:       Optional additional antenna gain (dBi).
                    To be used if not included in the pattern (ie when using
                    a normalized pattern).

  Returns:
    The antenna gains (in dB): either a scalar if hor_dirs is scalar,
    or an ndarray otherwise.
  """
  is_scalar = np.isscalar(hor_dirs)
  hor_dirs = np.atleast_1d(hor_dirs)

  bore_angle = hor_dirs - ant_azimuth
  bore_angle[bore_angle >= 360] -= 360
  bore_angle[bore_angle < 0] += 360
  idx0 = bore_angle.astype(np.int)
  alpha = bore_angle - idx0
  idx1 = idx0 + 1
  idx1[idx1 >= 360] -= 360
  gains = (1-alpha) * hor_pattern[idx0] + alpha * hor_pattern[idx1]
  gains += ant_gain

  if is_scalar: return gains[0]
  return gains


def GetStandardAntennaGains(hor_dirs, ant_azimuth=None, ant_beamwidth=None, ant_gain=0):
  """Computes the antenna gains from a standard antenna defined by beamwidth.

  See R2-SGN-20.
  This uses the standard 3GPP formula for pattern derivation from a given
  antenna 3dB cutoff beamwidth.
  Directions and azimuth are defined compared to the north in clockwise
  direction and shall be within [0..360] degrees.

  Inputs:
    hor_dirs:       Ray directions in horizontal plane (degrees).
                    Either a scalar or an iterable.
    ant_azimut:     Antenna azimuth (degrees).
    ant_beamwidth:  Antenna 3dB cutoff beamwidth (degrees).
                    If None, then antenna is isotropic (default).
    ant_gain:       Antenna gain (dBi).

  Returns:
    The CBSD antenna gains (in dB).
    Either a scalar if hor_dirs is scalar or an ndarray otherwise.
  """
  is_scalar = np.isscalar(hor_dirs)
  hor_dirs = np.atleast_1d(hor_dirs)

  if (ant_beamwidth is None or ant_azimuth is None or
      ant_beamwidth == 0 or ant_beamwidth == 360):
    gains = ant_gain * np.ones(hor_dirs.shape)
  else:
    bore_angle = hor_dirs - ant_azimuth
    bore_angle[bore_angle > 180] -= 360
    bore_angle[bore_angle < -180] += 360
    gains = -12 * (bore_angle / float(ant_beamwidth))**2
    gains[gains < -20] = -20.
    gains += ant_gain

  if is_scalar: return gains[0]
  return gains


def GetRadarNormalizedAntennaGains(hor_dirs,
                                   radar_azimuth,
                                   radar_beamwidth=3):
  """Computes the DPA radar normalized antenna gain.

  See R2-SGN-24.
  Directions and azimuth are defined compared to the north in clockwise
  direction and shall be within [0..360] degrees.
  Note that the DPA antenna gain is normalized to 0dBi at boresight:
  actual radar antenna gain is implicitely included in the target interference
  thresholds.

  Inputs:
    hor_dirs:       Ray directions in horizontal plane (degrees).
                    Either a scalar or an iterable.
    radar_azimuth:  The radar antenna azimuth (degrees).
    radar_beamwidth: The radar antenna beamwidth (degrees).

  Returns:
    The normalized antenna gains (in dB).
    Either a scalar if hor_dirs is scalar or an ndarray otherwise.

  """
  if radar_beamwidth == 360:
    return 0.
  is_scalar = np.isscalar(hor_dirs)
  hor_dirs = np.atleast_1d(hor_dirs)

  bore_angle = hor_dirs - radar_azimuth
  bore_angle[bore_angle > 180] -= 360
  bore_angle[bore_angle < -180] += 360
  bore_angle = np.abs(bore_angle)
  gains = -25 * np.ones(len(bore_angle))
  gains[bore_angle < radar_beamwidth / 2.] = 0

  if is_scalar: return gains[0]
  return gains


def GetFssAntennaGains(hor_dirs, ver_dirs,
                       fss_pointing_azimuth, fss_pointing_elevation,
                       fss_antenna_gain,
                       w1=0, w2=1.0):
  """Computes the FSS earth station antenna gain.

  See R2-SGN-21.
  Horizontal directions and azimuth are defined compared to the north in
  clockwise fashion and shall be within [0..360] degrees.
  Vertical directions are positive for above horizon, and negative below.

  Inputs:
    hor_dirs:               Ray directions in horizontal plane (degrees).
                            Either a scalar or an iterable.
    ver_dirs:               Ray directions in vertical plane (degrees).
                            Either a scalar or an iterable (same dimension
                            as hor_dirs)
    fss_pointing_azimuth:   FSS earth station azimuth angle (degrees).
    fss_pointing_elevation: FSS earth station vertical angle (degrees).
    fss_antenna_gain:       FSS earth station nominal antenna gain (dBi).
    w1, w2:                 Weights on the tangent and perpendicular
                            components. Optional: by default only use the
                            perpendicular component.
  Returns:
    The FSS gains on the incoming ray (in dB).
    Either a scalar if hor_dirs is scalar, or an ndarray otherwise.
  """
  is_scalar = np.isscalar(hor_dirs)
  hor_dirs = np.atleast_1d(np.radians(hor_dirs))
  ver_dirs = np.atleast_1d(np.radians(ver_dirs))
  fss_pointing_elevation = np.radians(fss_pointing_elevation)
  fss_pointing_azimuth = np.radians(fss_pointing_azimuth)

  # Compute the satellite antenna off-axis angle - see formula in R2-SGN-21, iii
  theta = 180/np.pi * np.arccos(
      np.cos(ver_dirs) * np.cos(fss_pointing_elevation)
      * np.cos(fss_pointing_azimuth - hor_dirs) +
      np.sin(ver_dirs) * np.sin(fss_pointing_elevation))

  gain_gso_t, gain_gso_p = _GetGsoGains(theta, fss_antenna_gain)
  gains = w1 * gain_gso_t + w2 * gain_gso_p

  if is_scalar: return gains[0]
  return gains


def _GetGsoGains(theta, nominal_gain):
  """Returns FSS earth station gains from the off-axis angle.

  GSO means the 'Geostationary Satellite Orbit'.

  Inputs:
    theta:        Off-axis angles (degrees), as a ndarray
    nominal_gain: Nominal antenna gain (dBi)
  Returns:
    a tuple of ndarray:
      gain_gso_t: Gains in the tangent plane of the GSO (dB).
      gain_gso_p: Gains in the perpendicular plane of the GSO (dB).
  """
  theta = np.abs(theta)
  gain_gso_t = -10 * np.ones(len(theta))
  gain_gso_p = gain_gso_t.copy()

  gain_gso_p[theta <= 3] = nominal_gain
  idx_3_to_48 = np.where((theta > 3) & (theta <= 48))[0]
  gain_gso_p[idx_3_to_48] = 32 - 25 * np.log10(theta[idx_3_to_48])

  gain_gso_t[theta <= 1.5] = nominal_gain
  idx_1_5_to_7 = np.where((theta > 1.5) & (theta <= 7))[0]
  gain_gso_t[idx_1_5_to_7] = 29 - 25 * np.log10(theta[idx_1_5_to_7])
  gain_gso_t[(theta > 7) & (theta <= 9.2)] = 8
  idx_9_2_to_48 = np.where((theta > 9.2) & (theta <= 48))[0]
  gain_gso_t[idx_9_2_to_48] = 32 - 25 * np.log10(theta[idx_9_2_to_48])

  return gain_gso_t, gain_gso_p
