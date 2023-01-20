"""
Robust calibration
"""
import numpy as np


def robviscal(vis_src, sigv_nr=0):
    """
    Robust calibration using SVD
    """
    # Compute rank one approx of source centered visibilities
    U_xx, S_xx, Vh_xx = np.linalg.svd(vis_src[0, 0, :, :], hermitian=True)
    U_xy, S_xy, Vh_xy = np.linalg.svd(vis_src[0, 1, :, :], hermitian=False)
    U_yx, S_yx, Vh_yx = np.linalg.svd(vis_src[1, 0, :, :], hermitian=False)
    U_yy, S_yy, Vh_yy = np.linalg.svd(vis_src[1, 1, :, :], hermitian=True)
    vis_src_sv_xx = S_xx[sigv_nr]*np.outer(U_xx[:, sigv_nr], Vh_xx[sigv_nr, :])
    vis_src_sv_xy = S_xy[sigv_nr]*np.outer(U_xy[:, sigv_nr], Vh_xy[sigv_nr, :])
    vis_src_sv_yx = S_yx[sigv_nr]*np.outer(U_yx[:, sigv_nr], Vh_yx[sigv_nr, :])
    vis_src_sv_yy = S_yy[sigv_nr]*np.outer(U_yy[:, sigv_nr], Vh_yy[sigv_nr, :])
    vis_src_sv = np.array([[vis_src_sv_xx, vis_src_sv_xy],
                           [vis_src_sv_yx, vis_src_sv_yy]])
    
    # Solve for gains
    gyrx = np.exp(0*0.55*2.0j*np.pi)
    vis1vecX = 1.00*np.sqrt(S_xx[sigv_nr])*U_xx[:,sigv_nr]
    vis1vecY = gyrx*np.sqrt(S_yy[sigv_nr])*U_yy[:,sigv_nr]
    #g_est_inv = np.array([np.conj(vis1vecX), np.conj(vis1vecY)])  # Just phase
    g_est_inv = np.array([1/vis1vecX, 1/vis1vecY])  # Amp & phase
    #g_est_inv = np.ones((2,vis1vecX.shape[0]))  # TEST no gain cal
    normfac = np.sqrt(np.sum(np.abs(g_est_inv)**2))
    g_est_inv = g_est_inv/normfac

    return g_est_inv
