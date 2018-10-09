import numpy as np
from SN_Rate import SN_Rate


class Generate_Sample:
    """ Generates a sample of parameters for simulation
    Input
    ---------
    sn_parameters: X1, Color, z, DayMax, ...
    (see yaml file input)
    cosmology parameters: H0, Om0
    Returns (call)
    ---------
    recordarray of parameters to simulate:
    X1,Color,z,DayMax
    """

    def __init__(self, sn_parameters, cosmo_parameters):

        self.params = sn_parameters
        self.sn_rate = SN_Rate(rate=self.params['z']['rate'],
                               H0=cosmo_parameters['H0'],
                               Om0=cosmo_parameters['Omega_m'])

        self.x1_color = self.Get_Dist(self.params['X1_Color']['rate'])

    def __call__(self, obs):
        """
        Input
        ---------
        recarray of observations

        Returns
        ---------
        recarray of sn parameters for simulation:
        z, X1 , Color ,  DayMax
        """

        # get duration of obs
        duration = np.max(obs['mjd'])-np.min(obs['mjd'])
        # get z range
        zmin = self.params['z']['min']
        zmax = self.params['z']['max']
        # get sn rate for this z range
        zz, rate, err_rate, nsn, err_nsn = self.sn_rate(
            zmin=zmin, zmax=zmax,
            duration=duration,
            survey_area = np.unique(obs['pixarea']))
        # get number of supernovae
        N_SN = int(np.cumsum(nsn)[-1])
        weight_z = np.cumsum(nsn)/np.sum(np.cumsum(nsn))
        dist_daymax = np.arange(np.min(obs['mjd']), np.max(obs['mjd']), 0.1)
        r = []
        for j in range(N_SN):
            z = self.Get_Val(self.params['z']['type'], zmin, zz, weight_z)
            x1_color = self.Get_Val(self.params['X1_Color']['type'],
                                    self.params['X1_Color']['min'],
                                    self.x1_color['low_z'][['X1', 'Color']],
                                    self.x1_color['low_z']['weight'])
            DayMax = self.Get_Val(self.params['DayMax']['type'],
                                  self.params['DayMax']['min'], dist_daymax,
                                  [1./len(dist_daymax)]*len(dist_daymax))
            r.append((z, x1_color[0], x1_color[1], DayMax))

        return np.rec.fromrecords(r, names=['z', 'X1', 'Color', 'DayMax'])

    def Get_Val(self, type, val, distrib, weight):
        """ Get values of a given parameter
        Input
        ---------
        type: random or not
        val : return value (if not random)
        distrib, weight: distrib and weight to get the parameter

        Returns
        ---------
        parameter value
        """

        if type == 'random':
            return np.random.choice(distrib, 1, p=weight)[0]
        else:
            return val

    def Get_Dist(self, rate):
        """ get (X1,C) distributions
        Input
        ---------
        rate: name of the X1_C distrib (JLA, ...)

        Returns
        ---------
        dict of (X1,C) rates
        keys: 'low_z' and 'high_z'
        val: recarray with X1,Color,weight_X1,weight_Color,weight
        """

        prefix = '../SN_Utils/input/Dist_X1_Color_'+rate+'_'
        suffix = '.txt'
        # names=['x1','c','weight_x1','weight_c','weight_tot']
        dtype = np.dtype([('X1', np.float), ('Color', np.float),
                          ('weight_X1', np.float), ('weight_Color', np.float),
                          ('weight', np.float)])
        x1_color = {}
        for val in ['low_z', 'high_z']:
            x1_color[val] = np.loadtxt(prefix+val+suffix, dtype=dtype)

        return x1_color

    def Plot_Parameters(self, gen_params):
        """ Plot the generated parameters
        (z,X1,Color,DayMax)
        """
        import pylab as plt

        fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(10, 9))

        todraw = ['z', 'X1', 'Color', 'DayMax']
        idx = dict(zip(todraw, [(0, 0), (0, 1), (1, 1), (1, 0)]))

        for name in todraw:
            i = idx[name][0]
            j = idx[name][1]
            ax[i][j].hist(gen_params[name], histtype='step')
            ax[i][j].set_xlabel(name)
            ax[i][j].set_ylabel('Number of entries')

        plt.show()
