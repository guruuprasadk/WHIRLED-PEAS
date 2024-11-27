# Copyright (c) 2014, Dignity Health
# 
#     The GPI core node library is licensed under
# either the BSD 3-clause or the LGPL v. 3.
# 
#     Under either license, the following additional term applies:
# 
#         NO CLINICAL USE.  THE SOFTWARE IS NOT INTENDED FOR COMMERCIAL
# PURPOSES AND SHOULD BE USED ONLY FOR NON-COMMERCIAL RESEARCH PURPOSES.  THE
# SOFTWARE MAY NOT IN ANY EVENT BE USED FOR ANY CLINICAL OR DIAGNOSTIC
# PURPOSES.  YOU ACKNOWLEDGE AND AGREE THAT THE SOFTWARE IS NOT INTENDED FOR
# USE IN ANY HIGH RISK OR STRICT LIABILITY ACTIVITY, INCLUDING BUT NOT LIMITED
# TO LIFE SUPPORT OR EMERGENCY MEDICAL OPERATIONS OR USES.  LICENSOR MAKES NO
# WARRANTY AND HAS NOR LIABILITY ARISING FROM ANY USE OF THE SOFTWARE IN ANY
# HIGH RISK OR STRICT LIABILITY ACTIVITIES.
# 
#     If you elect to license the GPI core node library under the LGPL the
# following applies:
# 
#         This file is part of the GPI core node library.
# 
#         The GPI core node library is free software: you can redistribute it
# and/or modify it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version. GPI core node library is distributed
# in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
# the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU Lesser General Public License for more details.
# 
#         You should have received a copy of the GNU Lesser General Public
# License along with the GPI core node library. If not, see
# <http://www.gnu.org/licenses/>.


# Author: Guru Krishnamoorthy
# Date: 2024Sep26

import gpi
import numpy as np


class ExternalNode(gpi.NodeAPI):
    '''Node summary goes here.
    '''

    def initUI(self):
        
        self.addWidget('PushButton', 'compute', val=False, toggle=True)
        
        self.addWidget('DoubleSpinBox', 'Global Delay(us)', val=0)
        self.addWidget('DoubleSpinBox', 'X Delay(us)', val=0)
        self.addWidget('DoubleSpinBox', 'Y Delay(us)', val=0)
        
        # IO Ports
        self.addInPort('params_in', 'PASS')
        
        self.addOutPort('crds_in_gre', 'NPYarray')
        self.addOutPort('sdc_in_gre', 'NPYarray')
        
        self.addOutPort('crds_out_gre', 'NPYarray')
        self.addOutPort('sdc_out_gre', 'NPYarray')
        
        self.addOutPort('crds_in_se', 'NPYarray')
        self.addOutPort('sdc_in_se', 'NPYarray')
        
        self.addOutPort('crds_out_se', 'NPYarray')
        self.addOutPort('sdc_out_se', 'NPYarray')
        
        
        self.addOutPort('grad_gre', 'NPYarray')
        self.addOutPort('grad_se', 'NPYarray')
        
        
    def compute(self):
        
        inparam = self.getData('params_in')
        if inparam is not None:
            
            psd_name = inparam['psd_iname']
            
            if self.getVal('compute'):
                if psd_name == 'SPLQSE':
                    [crds_in_se, crds_out_se, sdc_in_se, sdc_out_se, 
                    crds_in_gre, crds_out_gre, sdc_in_gre, sdc_out_gre,
                    grad_se, grad_gre] = gen_traj_splqse(inparam)
                    
                    global_delay = self.getVal('Global Delay(us)')
                    x_delay = self.getVal('X Delay(us)')
                    y_delay = self.getVal('Y Delay(us)')
                    
                    
                    dwell = inparam['sp_gen_dwell'] #usec
                    
                    
                    if crds_in_se is not None:
                        crds_in_se = add_delay(crds_in_se,  global_delay, x_delay, y_delay, dwell)
                        
                    if crds_out_se is not None:
                        crds_out_se = add_delay(crds_out_se,  global_delay, x_delay, y_delay, dwell)
                        
                    if crds_in_gre is not None:
                        crds_in_gre = add_delay(crds_in_gre,  global_delay, x_delay, y_delay, dwell)
                        
                    if crds_out_gre is not None:
                        crds_out_gre = add_delay(crds_out_gre,  global_delay, x_delay, y_delay, dwell)
                    
                    self.setData('crds_in_se', crds_in_se)
                    self.setData('crds_out_se', crds_out_se)
                    self.setData('sdc_in_se', sdc_in_se)
                    self.setData('sdc_out_se', sdc_out_se)
                    self.setData('crds_in_gre', crds_in_gre)
                    self.setData('crds_out_gre', crds_out_gre)
                    self.setData('sdc_in_gre', sdc_in_gre)
                    self.setData('sdc_out_gre', sdc_out_gre)
                    
                    self.setData('grad_se', grad_se)
                    self.setData('grad_gre', grad_gre)
                    
                if psd_name == 'SPIRALGRE3D':
                    [crds_in_gre, crds_out_gre, sdc_in_gre, sdc_out_gre,
                    grad_gre] = gen_traj_spgre3d(inparam)
                    
                    self.setData('crds_in_gre', crds_in_gre)
                    self.setData('crds_out_gre', crds_out_gre)
                    self.setData('sdc_in_gre', sdc_in_gre)
                    self.setData('sdc_out_gre', sdc_out_gre)
                    self.setData('grad_gre', grad_gre)
                    
        return 0
"""
Support functions
"""  

#%% # Traj Gen for SPGRE3D
def gen_traj_spgre3d(inparam):
        # Initialize parameters for SE trajectory generation
    params = {
        'fov': inparam['sp_gen_fov'],
        'res': inparam['sp_gen_res'],
        'm_slew': inparam['sp_gen_smax'],
        'm_grad': inparam['sp_gen_gmax'],
        'm_omega': inparam['sp_gen_omax'],
        'grast': inparam['sp_gen_rast'],
        'dwell': inparam['sp_gen_dwell'],
        'gamma': inparam['sp_gen_gamma'],
        'add_outer_ring': inparam['sp_do_outring'],
        'nyq_arms': inparam['sp_nyq_arms'],
        'dixon_acqs': 1, # later
        'n_intra_tr_echoes': inparam['sp_gen_intratr_necho'],
        'arm_order' : inparam['sp_gen_arm_order']
    }
    se_direction = inparam['sp_gen_direction']
    
    # Generate SE trajectory based on the direction
    crds_in_gre = None
    crds_out_gre = None
    sdc_in_gre = None
    sdc_out_gre = None
    grad_gre = None
    if se_direction == 'out':
        params.update({
            'req_tau': inparam['sp_gen_tau'],
            'spiral_direction': 'Spiral Out'
        })
        grad_gre, ksp_gre, sdc_gre, params_gre_out = wpgen_single_tau(params)
        [crds_out_gre, sdc_out_gre] = gen_rotated_traj(ksp_gre, sdc_gre, params)
    elif se_direction == 'in-out':
        params.update({
            'req_tau': inparam['sp_gen_tau'],
            'spiral_direction': 'Spiral In-Out'
        })
        grad_gre, ksp_gre, sdc_gre, params_gre_out = wpgen_single_tau(params)
        [crds_inout_gre, sdc_inout_gre] = gen_rotated_traj(ksp_gre, sdc_gre, params)
        
        #Split in and out
        gre_in_pts = crds_inout_gre.shape[-2] // 2
        crds_in_gre = crds_inout_gre[..., 0:gre_in_pts, :]
        crds_out_gre = crds_inout_gre[..., gre_in_pts:, :]
        
        sdc_in_gre = sdc_inout_gre[...,0:gre_in_pts]
        sdc_out_gre = sdc_inout_gre[...,gre_in_pts:]
        
    # Return the generated SE and GRE trajectories
    return (
        crds_in_gre, crds_out_gre, sdc_in_gre, sdc_out_gre, grad_gre
    )

#%% # WPGEN single Tau
def wpgen_single_tau(params):
    import guru.sample.cppHelper as cpphelper
    
    fov = params['fov']
    res = params['res']
    req_tau = params['req_tau']
    m_slew = params['m_slew']
    m_grad = params['m_grad']
    m_omega = params['m_omega']
    grast = params['grast']
    dwell = params['dwell']
    gamma = params['gamma']
    add_outer_ring = params['add_outer_ring']
    n_intra_tr_echoes = params['n_intra_tr_echoes']
    spiral_direct = {
                    'Spiral Out': cpphelper.SpiralDirect.SPIRAL_OUT,
                    'Spiral In': cpphelper.SpiralDirect.SPIRAL_IN,
                    'Spiral In-Out': cpphelper.SpiralDirect.SPIRAL_INOUT,
                    'Spiral Out-In': cpphelper.SpiralDirect.SPIRAL_OUTIN
                }.get(params['spiral_direction'], cpphelper.SpiralDirect.SPIRAL_OUT)
    
    wpgen = cpphelper.WPGen(fov, res, req_tau, m_slew, m_grad, m_omega, grast, dwell, gamma,
                                    add_outer_ring)
    
    wpcompose = cpphelper.WPCompose(wpgen, n_intra_tr_echoes, spiral_direct)
    
    g_out = wpcompose.composeGradients()
    [ksp, sdc] = wpcompose.composeKSPnSDC()
    
    out_params = {
        'time_echo1': wpcompose.getTimeEcho1(),
        'delta_te': wpcompose.getDeltaTE(),
        't_arc': wpgen.GetArcDur(),
        't_omega': wpgen.GetFreqDur(),
        't_slew': wpgen.GetSlewDur(),
        't_grad': wpgen.GetGradDur(),
        'tau_total': wpgen.GetActualTau(),
        't_ramp': wpgen.GetRampDur(),
        't_ring': wpgen.GetOutRingDur(),
        'omega_max': wpgen.GetOmegaMax(),
        'slew_max': wpgen.GetSlewMax(),
        'grad_max': wpcompose.getGradMax(),
        'snr_factor': wpgen.GetSnrFactor(),
        'ksp_pts': wpgen.GetKspacePts(),
        'nyq_arms': wpcompose.getNyqArms(),
        'grd_mtx': wpgen.GetGridMtxSize()
    }
    
    return g_out, ksp, sdc, out_params

#%% # WPGEN dual Tau
def wpgen_dual_tau(params):
    import guru.sample.cppHelper as cpphelper
    
    fov = params['fov']
    res = params['res']
    req_tau_in = params['req_tau_in']
    req_tau_out = params['req_tau_out']
    m_slew = params['m_slew']
    m_grad = params['m_grad']
    m_omega = params['m_omega']
    grast = params['grast']
    dwell = params['dwell']
    gamma = params['gamma']
    add_outer_ring = params['add_outer_ring']
    nav_pts = int(params['nav_pts'] / int(dwell * 1000))
    
    wpgen_in = cpphelper.WPGen(fov, res, req_tau_in, m_slew, m_grad, m_omega, grast, dwell, gamma,
                                    add_outer_ring)
    
    wpgen_out = cpphelper.WPGen(fov, res, req_tau_out, m_slew, m_grad, m_omega, grast, dwell, gamma,
                                    add_outer_ring)
    
    wpcompose = cpphelper.WPCompose(wpgen_in, wpgen_out, nav_pts)
    
    g_out = wpcompose.composeGradients()
    [ksp, sdc] = wpcompose.composeKSPnSDC()
    
    def extract_params(wpgen, wpcompose):
        return {
            'time_echo1': wpcompose.getTimeEcho1(),
            'delta_te': wpcompose.getDeltaTE(),
            'grad_max': wpgen.GetGradMax(),
            'nyq_arms': wpgen.GetNyqNumberArms(),
            't_arc': wpgen.GetArcDur(),
            't_omega': wpgen.GetFreqDur(),
            't_slew': wpgen.GetSlewDur(),
            't_grad': wpgen.GetGradDur(),
            'tau_total': wpgen.GetActualTau(),
            't_ramp': wpgen.GetRampDur(),
            't_ring': wpgen.GetOutRingDur(),
            'omega_max': wpgen.GetOmegaMax(),
            'slew_max': wpgen.GetSlewMax(),
            'snr_factor': wpgen.GetSnrFactor(),
            'ksp_pts': wpgen.GetKspacePts(),
            'grd_mtx': wpgen.GetGridMtxSize()
        }

    params_in = extract_params(wpgen_in, wpcompose)
    params_out = extract_params(wpgen_out, wpcompose)
    
    return g_out, ksp, sdc, params_in, params_out

#%% # Traj Gen for SPLQSE
def gen_traj_splqse(inparam):
    # Initialize parameters for SE trajectory generation
    params = {
        'fov': inparam['sp_gen_fov'],
        'res': inparam['sp_gen_res'],
        'm_slew': inparam['sp_gen_smax'],
        'm_grad': inparam['sp_gen_gmax'],
        'm_omega': inparam['sp_gen_omax'],
        'grast': inparam['sp_gen_rast'],
        'dwell': inparam['sp_gen_dwell'],
        'gamma': inparam['sp_gen_gamma'],
        'add_outer_ring': inparam['sp_do_outring'],
        'nyq_arms': inparam['sp_arms'],
        'dixon_acqs': 1, # later
        'n_intra_tr_echoes': 1,
        'arm_order' : inparam['sp_gen_arm_order']
    }
    se_direction = inparam['sp_gen_se_direction']

    # Generate SE trajectory based on the direction
    crds_in_se = None
    crds_out_se = None
    sdc_in_se = None
    sdc_out_se = None
    grad_se = None
    if se_direction == 'out':
        params.update({
            'req_tau': inparam['sp_gen_se_tau_out'],
            'spiral_direction': 'Spiral Out'
        })
        grad_se, ksp_se, sdc_se, params_se_out = wpgen_single_tau(params)
        [crds_out_se, sdc_out_se] = gen_rotated_traj(ksp_se, sdc_se, params)
    elif se_direction == 'in-out':
        params.update({
            'req_tau_in': inparam['sp_gen_se_tau_in'],
            'req_tau_out': inparam['sp_gen_se_tau_out'],
            'nav_pts': inparam['sp_gen_inout_nav_pts']
        })
        grad_se, ksp_se, sdc_se, params_se_in, params_se_out  = wpgen_dual_tau(params)
        [crds_inout_se, sdc_inout_se] = gen_rotated_traj(ksp_se, sdc_se, params)
        
        #Split in and out
        se_in_pts = params_se_in['ksp_pts']
        crds_in_se = crds_inout_se[..., 0:se_in_pts, :]
        crds_out_se = crds_inout_se[..., se_in_pts+params['nav_pts']:, :]
        
        sdc_in_se = sdc_inout_se[...,0:se_in_pts]
        sdc_out_se = sdc_inout_se[...,se_in_pts+params['nav_pts']:]

    # Initialize GRE trajectory as None
    crds_in_gre = None
    crds_out_gre = None
    sdc_in_gre = None
    sdc_out_gre = None
    gre_sp_direction = inparam['sp_gen_add_gre']
    grad_gre = None
    
    # Generate GRE trajectory if specified
    if(gre_sp_direction != 'off'):
        if(gre_sp_direction == 'out'):
            params.update({
                'req_tau': inparam['sp_gen_gre_tau_out'],
                'spiral_direction': 'Spiral Out'
            })
            grad_gre, ksp_gre, sdc_gre, params_gre_out = wpgen_single_tau(params)
            [crds_out_gre, sdc_out_gre] = gen_rotated_traj(ksp_gre, sdc_gre, params)
        elif(gre_sp_direction == 'in'):
            params.update({
                'req_tau': inparam['sp_gen_gre_tau_in'],
                'spiral_direction': 'Spiral In'
            })
            grad_gre, ksp_gre, sdc_gre, params_gre_in = wpgen_single_tau(params)
            [crds_in_gre, sdc_in_gre] = gen_rotated_traj(ksp_gre, sdc_gre, params)
        elif(gre_sp_direction == 'in-out'):
            params.update({
                'req_tau_in': inparam['sp_gen_gre_tau_in'],
                'req_tau_out': inparam['sp_gen_gre_tau_out'],
                'nav_pts': 0
            })
            grad_gre, ksp_gre, sdc_gre, params_gre_in, params_gre_out = wpgen_dual_tau(params)
            [crds_inout_gre, sdc_inout_gre] = gen_rotated_traj(ksp_gre, sdc_gre, params)
            gre_in_pts = params_gre_in['ksp_pts']

            #Split in and out
            crds_in_gre = crds_inout_gre[...,0:gre_in_pts,:]
            crds_out_gre = crds_inout_gre[...,gre_in_pts::,:]
            
            sdc_in_gre = sdc_inout_gre[...,0:gre_in_pts]
            sdc_out_gre = sdc_inout_gre[...,gre_in_pts::]
        elif(gre_sp_direction == 'out-in'):
            params.update({
                'req_tau': inparam['sp_gen_gre_tau_in'],
                'spiral_direction': 'Spiral Out-In'
            })
            grad_gre, ksp_gre, sdc_gre, params_gre_outin = wpgen_single_tau(params)
            [crds_inout_gre, sdc_inout_gre] = gen_rotated_traj(ksp_gre, sdc_gre, params)
            gre_pts = params_gre_outin['ksp_pts']

            #Split in and out
            crds_out_gre = crds_inout_gre[...,0:gre_pts,:]
            crds_in_gre = crds_inout_gre[...,gre_pts::,:]
            
            sdc_out_gre = sdc_inout_gre[...,0:gre_pts]
            sdc_in_gre = sdc_inout_gre[...,gre_pts::]
        
    # Return the generated SE and GRE trajectories
    return (
        crds_in_se, crds_out_se, sdc_in_se, sdc_out_se, 
        crds_in_gre, crds_out_gre, sdc_in_gre, sdc_out_gre, grad_se, grad_gre
    )

#%% # Traj Gen for SPLQSE
def gen_rotated_traj(ksp,sdc, params):
    import guru.sample.cppHelper as cpphelper
    
    arm_ord = {
                'Linear': cpphelper.OrderType.LINEAR,
                'Skip': cpphelper.OrderType.SKIP,
                'TwoWay': cpphelper.OrderType.TWO_WAY,
                'Mixed': cpphelper.OrderType.MIXED,
                'Golden': cpphelper.OrderType.GOLDEN
            }.get(params['arm_order'], cpphelper.OrderType.LINEAR)
    
    #angular_coverage = np.pi if spiral_direction == 'Spiral In-Out' else 2*np.pi
    angular_coverage = 2*np.pi
    sporder = cpphelper.SpiralArmOrder(params['nyq_arms'],  params['dixon_acqs'], angular_coverage,arm_ord,  False, False)
    arm_angles = sporder.composeAngles()

    crd_full = rotate_trajectory(ksp, arm_angles, params['nyq_arms'], params['n_intra_tr_echoes'], params['dixon_acqs'])
    sdc_full = np.broadcast_to(sdc, (params['dixon_acqs'], params['n_intra_tr_echoes'], params['nyq_arms'], sdc.shape[-1]))
    sdc_full = np.squeeze(sdc_full)
    
    return crd_full, sdc_full

#%% # Rotate k-space trajectory
def rotate_trajectory(arm_to_rotate, arm_angles, total_arms, total_echo, total_acq):
    sample_pts = arm_to_rotate.shape[-2]
    coords_full = np.zeros([total_acq, total_echo, total_arms, sample_pts, 2], dtype=np.double)
    
    for arm in range(total_arms):
        for echo in range(total_echo):
            for acq in range(total_acq):
                curr_angle = arm_angles[acq, arm]
                sb = np.sin(curr_angle)
                cb = np.cos(curr_angle)
                coords_full[acq, echo, arm, :, 0] = (cb * arm_to_rotate[echo, :, 0]) - (sb * arm_to_rotate[echo, :, 1])
                coords_full[acq, echo, arm, :, 1] = (cb * arm_to_rotate[echo, :, 1]) + (sb * arm_to_rotate[echo, :, 0])
    
    coords_full = np.squeeze(coords_full)
    return coords_full

    
  
#%% # Delay k-space trajectory  
def add_delay(crds_full,  global_delay, x_delay, y_delay,dwell):
    from scipy.interpolate import interp1d

    dwell *= 1e-3 #msec to usec
    nsample = crds_full.shape[-2]
    total_arms =  crds_full.shape[-3]
    tau = dwell * nsample
    crds_delayed = np.zeros_like(crds_full)
    
    
    t_samp = np.arange(0, tau, dwell) 
        
    
    shifted_time_x = t_samp + (x_delay * 1e-6)
    shifted_time_y = t_samp + (y_delay * 1e-6)
    
    shifted_time_x = shifted_time_x + (global_delay * 1e-6)
    shifted_time_y = shifted_time_y + (global_delay * 1e-6)
    
    for j in range(total_arms):
        
        kx = crds_full[j,:,0]
        ky = crds_full[j,:,1]

        kx_interpolator = interp1d(t_samp, kx, kind='quadratic', bounds_error=False, fill_value=(0, 0))
        ky_interpolator = interp1d(t_samp, ky, kind='quadratic', bounds_error=False, fill_value=(0, 0))
        
        crds_delayed[j,:,0] = kx_interpolator(shifted_time_x)
        crds_delayed[j,:,1] = ky_interpolator(shifted_time_y)

    return crds_delayed  