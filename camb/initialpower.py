# Initial power spectrum parameters

from .baseconfig import CAMB_Structure, CAMBError, camblib
from ctypes import c_int, c_double, c_bool, POINTER, byref
import numpy as np

nnmax = 5

tensor_param_indeptilt = 1
tensor_param_rpivot = 2
tensor_param_AT = 3


class InitialPowerParams(CAMB_Structure):
    """
    Object to store parameters for the primordial power spectrum.
    Many of the internal variables are arrays, to allow calculating more than one power spectrum at one. Higher-level functions in the
    CAMB python wrapper assume only one is calculated.

    """
    _fields_ = [
        ("tensor_parameterization", c_int),
        ("nn", c_int),
        ("an", c_double * nnmax),
        ("n_run", c_double * nnmax),
        ("n_runrun", c_double * nnmax),
        ("ant", c_double * nnmax),
        ("nt_run", c_double * nnmax),
        ("rat", c_double * nnmax),
        ("k_0_scalar", c_double),
        ("k_0_tensor", c_double),
        ("ScalarPowerAmp", c_double * nnmax),
        ("TensorPowerAmp", c_double * nnmax),
        ("log10f",c_double * nnmax),
        ("deltans",c_double * nnmax),
        ("deltaPhi",c_double * nnmax),
        ("paxionf",c_double * nnmax),
        ("p_axion", c_double),
        ("phiEND", c_double),
        ("phiZERO", c_double),
        ("NEfoldStar", c_double)
    ]

    def set_params(self, As=2e-9, ns=0.96, nrun=0, nrunrun=0, r=0, nt=None, ntrun=0,
                   pivot_scalar=0.05, pivot_tensor=0.05, log10_f = 0, delta_ns = 0, delta_Phi = 0,
                   paxion_f = 0, paxion = 4./3., phi_END  = 0.59, phi_ZERO = 12.38, NEfold_Star = 57.5,
                   parameterization=tensor_param_rpivot):
        """
        Set parameters using standard power law parameterization. If nt=None, uses inflation consistency relation.

        :param As: comoving curvature power at k=piveo_scalar
        :param ns: scalar spectral index
        :param nrun: running of scalar spectral index d n_s/d log k
        :param nrunrun: running of running of spectral index
        :param r: tensor to scalar ratio at pivot
        :param nt: tensor spectral index. If None, set using inflation consistency
        :param ntrun: running of tensor spectral index
        :param pivot_scalar: pivot scale for scalar spectrum
        :param pivot_tensor:  pivot scale for tensor spectrum
        :param log10_f:  frequency
        :param delta_ns:  amplitude axion monodromy power spectrum
        :param delta_Phi:  phase
        :param paxion_f:  powerlaw index axion power spectrum
        :param paxion: see notes
        :param phi_END: see notes
        :param phi_ZERO: see notes
        :param NEfold_Star: see notes
        :param parameterization: See CAMB notes. One of
            - tensor_param_indeptilt = 1
            - tensor_param_rpivot = 2
            - tensor_param_AT = 3
        :return: self
        """

        if not parameterization in [tensor_param_rpivot, tensor_param_indeptilt]:
            raise CAMBError('Initial power parameterization not supported here')
        self.tensor_parameterization = parameterization
        self.ScalarPowerAmp[0] = As
        self.nn = 1  # number of different power spectra
        self.an[0] = ns
        self.n_run[0] = nrun
        self.n_runrun[0] = nrunrun
        self.log10f[0] = log10_f
        self.deltans[0] = delta_ns
        self.deltaPhi[0] = delta_Phi
        self.paxionf[0] = paxion_f
        if nt is None:
            # set from inflationary consistency
            if ntrun: raise CAMBError('ntrun set but using inflation consistency (nt=None)')
            if tensor_param_rpivot != tensor_param_rpivot:
                raise CAMBError('tensor parameterization not tensor_param_rpivot with inflation consistency')
            self.ant[0] = - r / 8.0 * (2.0 - ns - r / 8.0)
            self.nt_run[0] = r / 8.0 * (r / 8.0 + ns - 1)
        else:
            self.ant[0] = nt
            self.nt_run[0] = ntrun
        self.rat[0] = r
        self.k_0_scalar = pivot_scalar
        self.k_0_tensor = pivot_tensor
        self.p_axion = paxion
        self.phiEND = phi_END
        self.phiZERO = phi_ZERO
        self.NEfoldStar = NEfold_Star
        return self

    def has_tensors(self):
        """
        Do these settings have non-zero tensors?

        :return: True of non-zero tensor amplitude
        """
        return self.rat[0]
