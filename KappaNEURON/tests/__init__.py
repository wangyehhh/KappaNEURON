import os
import pkgutil
import tempfile
import shutil
import KappaNEURON
import unittest
import neuron
from neuron import *
from neuron import rxd
import numpy as np
import matplotlib.pyplot as plt
import re
import platform
import glob
from neuron.rxd.generalizedReaction import molecules_per_mM_um3

def compile_modfiles(dirpath='.'):
    cwd = os.getcwd()
    print("CWD: " + cwd)
    os.chdir(dirpath)
    print("CWD: " + os.getcwd())
    ## Find the version of nrnivmodl corresponding to this python module
    neuron_root = os.path.realpath(os.path.join(os.path.realpath(pkgutil.get_loader('neuron').filename), '../../../'))
    print("neuron_root: " + neuron_root)
    ## os.chdir(dirpath)
    ## Use it to compile files in the test directory
    ## os.system(os.path.join(neuron_root, 'x86_64/bin/nrnivmodl test'))
    if platform.system() == 'Windows':
        cmd = 'sh ' + os.path.realpath(os.path.join(neuron_root, 'bin/mknrndll'))
    else:
        cmd = os.path.realpath(os.path.join(neuron_root, 'x86_64/bin/nrnivmodl'))
    print("Compiling: " + cmd + "...")
    os.system(cmd)
    print('...compiled')
    os.chdir(cwd)
    print("CWD: " + cwd)

class TestCaAccumulation(unittest.TestCase):
    ## Compile and load mechanisms - note that the version of
    ## nrnivmodl must match NEURON dll - there is no way of checking
    ## this at present
    print(os.getcwd())
    if platform.system() != 'Windows':
        dirpath = tempfile.mkdtemp()
        for file in glob.glob(os.path.join(pkgutil.get_loader("KappaNEURON").filename, 'tests', '*.mod')):
            shutil.copy(file, dirpath)
    
        compile_modfiles(dirpath)
        ## neuron.load_mechanisms(dirpath)
        neuron.load_mechanisms(dirpath)
        ## We can't delete the mechanisms
        ## shutil.rmtree(dirpath)
    
        ## Whether to plot
    plot = True

    ## We can't put this stuff in setUp(), since problems arise if we
    ## try to redefine sections, regions and species. This is because
    ## of class variables in Species.
    ## Set up a kappa simulation of one compartment
    sk = h.Section()
    sk.insert("capulse")
    sk.L=0.2
    sk.diam=0.5
    r = rxd.Region([sk], nrn_region='i')
    ca = rxd.Species(r, name='ca', charge=2, initial=0.00005)

    ## Set up a mod simulation of one compartment
    sm = h.Section()
    sm.insert("capulse")
    sm.L=0.2
    sm.diam=0.5
    
    ## Time of pulse
    t0 = 1.0
    t1 = 1.1
    tstop = t1

    ## Conversions
    NA =  6.02214129e23     # Avogadro's number

    gbar = 0.001
    cm = sk(0.5).cm

    k1 = 0                      # Pump parameter
    k2 = False                  # Optional parameter
    ## P = False              # Pump species
    P = None
    P0 = 0
    
    tol = 0.01
    KappaNEURON.verbose = False
    mechanism = None

    def assertEqualWithinTol(self, a, b, tol=None):
        if tol == None:
            tol = self.tol
        self.assertAlmostEqual(a, b, delta=tol*a)

    def setUp(self):
        ## This runs at the start of every test...
        self.caitonum = molecules_per_mM_um3*np.pi*(self.sk.diam**2)/4*self.sk.L
        h.dt = 0.025/4
        print ("")
        print ("======================================================================")
        print ("In method", self._testMethodName)
        print ("----------------------------------------------------------------------")

    def get_mode(self, sec):
        ## Determine if section contains mod pump or kappa pump
        mode = 'kappa'
        for mech in sec(0.5):
            if re.search('caPump', mech.name()):
                mode = 'mod'
        return(mode)

    def injectCalcium(self, ghk=0, mechanism='caPump1'):
        print(pkgutil.get_loader('KappaNEURON').filename)
        print(os.path.realpath(KappaNEURON.__file__))
        ## Insert calcium pump into mod section
        self.sm.insert(mechanism)
        self.mechanism = mechanism

        ## Insert calcium pump into kappa section
        if mechanism == 'caPump1':
            print(KappaNEURON.__file__)            
            self.kappa = KappaNEURON.Kappa(membrane_species=[self.ca], kappa_file=os.path.dirname(KappaNEURON.__file__) + "/tests/" + mechanism + ".ka", regions=self.r)
        if mechanism == 'caPump2':
            self.P  = rxd.Species(self.r, name='P', charge=0, initial=self.P0)
            self.kappa = KappaNEURON.Kappa(membrane_species=[self.ca], species=[self.P], kappa_file=os.path.dirname(KappaNEURON.__file__) + "/tests/" + mechanism + ".ka", regions=self.r)
            self.kappa.setVariable('vol', self.sk.L*(self.sk.diam**2)/4*np.pi)
            self.kappa.setVariable('k2', self.k2)
            setattr(self.sm(0.5), 'k2_' + mechanism, self.k2)
            self.sm(0.5).P0_caPump2 = self.P0

        ## Set variables
        self.kappa.setVariable('k1', self.k1)
        setattr(self.sm(0.5), 'k1_' + mechanism, self.k1)

        ## Set up recordings
        self.rec_t = h.Vector()
        self.rec_t.record(h._ref_t)
        self.rec_cai = []
        self.rec_v = []
        for sec in h.allsec():
            self.rec_cai.append(h.Vector())
            self.rec_cai[-1].record(sec(0.5)._ref_cai)
            self.rec_v.append(h.Vector())
            self.rec_v[-1].record(sec(0.5)._ref_v)

        # self.sm(0.5).P0_caPump2 = self.P0
        if mechanism == 'caPump2':
            self.rec_Pi = []
            for sec in h.allsec():
                self.rec_Pi.append(h.Vector())
                print (self.get_mode(sec))
                if self.get_mode(sec) == 'mod':
                    self.rec_Pi[-1].record(sec(0.5)._ref_P_caPump2)
                else:
                    self.rec_Pi[-1].record(sec(0.5)._ref_Pi)
        
        KappaNEURON.progress = False

        self.assertIsInstance(self.sk, nrn.Section)

        for sec in h.allsec():
            ## This forces eca to be a constant, rather than being
            ## computed from Nernst equation at every time step
            ## See http://www.neuron.yale.edu/neuron/static/new_doc/modelspec/programmatic/ions.html?highlight=ion_style#ion_style
            if (ghk == 1):
                sec.fghk_capulse = 1
                h.ion_style("ca_ion", 3, 2, 1, 1, 1, sec=sec)
            else:
                h.ion_style("ca_ion", 3, 1, 0, 0, 1, sec=sec)
                
            for seg in sec:
                seg.t0_capulse = self.t0
                seg.t1_capulse = self.t1

                seg.gbar_capulse = self.gbar

        ## Initialise simulation
        ## neuron.init() does not take the v_init argument, so use:
        print("Before neuron.h.finitialize(-65.0) run")
        print("kappa section: ECa = %f; cai = %f; cao = %f; V = %f" % (self.sk(0.5).eca, self.sk(0.5).cai, self.sk(0.5).cao, self.sk(0.5).v))
        print("mod section:   ECa = %f; cai = %f; cao = %f; V = %f" % (self.sm(0.5).eca, self.sm(0.5).cai, self.sm(0.5).cao, self.sm(0.5).v))
        neuron.h.finitialize(-65.0)
        print("neuron.h.finitialize(-65.0) has been run")
        eca0 = self.sk(0.5).eca
        print("kappa section: ECa = %f; cai = %f; cao = %f; V = %f" % (self.sk(0.5).eca, self.sk(0.5).cai, self.sk(0.5).cao, self.sk(0.5).v))
        print("mod section:   ECa = %f; cai = %f; cao = %f; V = %f" % (self.sm(0.5).eca, self.sm(0.5).cai, self.sm(0.5).cao, self.sm(0.5).v))
        self.v0 = self.sk(0.5).v
        self.cai0 = self.sk(0.5).cai
        self.assertEqual(h.t, 0.0)

        ## Initialisation ensures that the concentrations is converted
        ## to a number of ions
        nions = round(0.00005*self.caitonum)
        self.assertEqual(self.sk(0.5).cai, nions/self.caitonum)

        ## FIXME: Should we need to run this *after* init()? Maybe
        ## this is an issue with the rxd.species code rather than
        ## ours, since ion_style() is called within the init sequence
        ## there.
        if (ghk == 0):
            h.ion_style("ca_ion", 3, 1, 0, 0, 1, sec=self.sk)

        ## Run
        run(self.tstop)
        print("After run()")
        print("kappa section: ECa = %f; cai = %f; cao = %f" % (self.sk(0.5).eca, self.sk(0.5).cai, self.sk(0.5).cao))
        print("mod section:   ECa = %f; cai = %f; cao = %f" % (self.sm(0.5).eca, self.sm(0.5).cai, self.sm(0.5).cao))

        self.assertAlmostEqual(h.t, self.tstop)
        self.assertGreater(self.sk(0.5).cai, 0.0)
        if (ghk == 0):
            self.assertEqual(eca0, self.sk(0.5).eca)

    def get_diffv_diffca(self, i):
        times = np.array(self.rec_t)
        stim_inds = np.where((times > self.t0) & (times < self.t1))
        ## Check that during the stimulus, every voltage increment
        ## is proportional to the calcium increment in the
        ## *preceeding* timestep.  At the end of the pulse and the
        ## beginning of the pulse this is not true, because the voltage increment 
        diffv = np.diff(np.array(self.rec_v[i])[stim_inds])
        diffca = np.diff(np.array(self.rec_cai[i])[stim_inds])
        return(diffv, diffca)

    def do_plot(self):
        plt.subplots_adjust(left=0.25)
        nrow = 3
        if (self.P):
            nrow = 4
        fig, ax = plt.subplots(nrows=nrow, ncols=1, figsize=(2.25*2, 2.5*2))
        fig.suptitle(self._testMethodName)
        i = 0
        labels=['kappa', 'mod']
        for sec in h.allsec():
            (diffv, diffca) = self.get_diffv_diffca(i)
            ax[0].plot(self.rec_t, self.rec_v[i], color='br'[i], label=labels[i])
            (Deltav_theo, Deltaca_theo) = self.get_Deltav_Deltaca_theo(sec, np.array(self.rec_t), False)
            ax[1].plot(self.rec_t, self.rec_cai[i], color='br'[i])
            ax[2].plot(diffv[1:len(diffv)-1], self.caitonum*diffca[0:len(diffv)-2], 'o', color='br'[i])
            if (self.P0 > 0):
                ax[3].plot(self.rec_t, self.rec_Pi[i], color='br'[i])
                ax[3].set_xlabel("time (ms)")
                ax[3].set_ylabel("P (mM)")
            fig.show()        
            i = i + 1

        ax[0].set_xlabel("time (ms)")
        ax[0].set_ylabel("V (mV)")
        ax[0].plot(self.rec_t, Deltav_theo + self.v0, color='g', label='theo')
        ax[0].legend(loc='upper left')
        ax[1].set_xlabel("time (ms)")
        ax[1].set_ylabel("Ca (mM)")
        ax[2].set_xlabel("Delta V")
        ax[2].set_ylabel("Delta #Ca ions")
        try:
            os.mkdir('test_figs')
        except:
            print('dir exists')
            
        plt.savefig(os.path.join('test_figs', self._testMethodName))

            
    def get_Deltav_Deltaca_theo(self, sec, t, output=True, verbose=False):
        eca = sec(0.5).eca
        volbyarea = sec.diam/4
        vtocai = self.cm/(1E-1*2*h.FARADAY*volbyarea)
        ## Print some variables
        v1 = sec(0.5).v
        if verbose:
            print("cca=%f, t0=%f, t1=%f, gbar=%f, cm=%f, v0=%f, v1=%f" % (sec(0.5).eca, self.t0, self.t1, self.gbar, self.cm, self.v0, v1))
        if verbose:
            print('Theoretical voltage and Ca difference:')
        tau = self.cm/(1000.0*self.gbar + self.k1*self.cm)
        vinf = (1000.0*self.gbar*eca + self.cm*self.k1*self.v0)/(1000.0*self.gbar + self.cm*self.k1)
        Deltav_theo = (vinf - self.v0)*(1-np.exp(-(t - self.t0)/tau))
        Deltav_theo1 = (vinf - self.v0)*(1-np.exp(-(self.t1 - self.t0)/tau))
        if isinstance(Deltav_theo, numpy.ndarray):
            Deltav_theo[np.where(t < self.t0)] = 0.0
        Deltaca_theo = Deltav_theo*vtocai
        Deltaca_theo1 = Deltav_theo1*vtocai
        if verbose:
            print("eca= %f; tau=%f; vinf=%f; Deltav_theo=%f; Deltaca_theo=%f" % (eca, tau, vinf, Deltav_theo, Deltaca_theo))
        if output:
            print('Theoretical voltage and Ca changes between start and end of injection:')
            print("  v(t1) -   v(t0) = %10.6f - %10.6f = %10.6f" % (self.v0 +   Deltav_theo1,  self.v0,   Deltav_theo1))
            print("cai(t1) - cai(t0) = %10.6f - %10.6f = %10.6f" % (self.cai0 + Deltaca_theo1, self.cai0, Deltaca_theo1))
            print()
        
        return(Deltav_theo, Deltaca_theo)

    def get_Deltav_Deltaca(self, sec, i):
        v1 = sec(0.5).v
        Deltav = v1 - self.v0
        Deltaca = sec(0.5).cai - self.rec_cai[i][0]
        print('Simulated voltage and Ca changes between start and end of injection:')
        print("  v(t1) -   v(t0) = %10.6f - %10.6f = %10.6f" % (v1, self.v0, Deltav))
        print("cai(t1) - cai(t0) = %10.6f - %10.6f = %10.6f" % (sec(0.5).cai, self.rec_cai[i][0], Deltaca))
        print()
        
        return(Deltav, Deltaca)

    def get_stats(self):
        Deltav       = {}
        Deltaca      = {}
        Deltav_theo  = {}
        Deltaca_theo = {}
        volbyarea    = {}
        vtocai       = {}
        diffv        = {}
        diffca       = {}
        
        i = 0
        for sec in h.allsec():
            ## Determine if section contains mod pump or kappa pump
            mode = self.get_mode(sec)
            print("----------------------------------------------------------------------")
            print("Section " + sec.name() + " contains " +  mode + " pump")
            print("----------------------------------------------------------------------")
            ## Print some variables
            v1 = sec(0.5).v
            print("t0=%f, t1=%f, gbar=%f, cm=%f" % (self.t0, self.t1, self.gbar, self.cm))
            print("Eca=%f, v0=%f, v1=%f" % (sec(0.5).eca, self.v0, v1))
            print
            (Deltav_theo[mode],  Deltaca_theo[mode])  = self.get_Deltav_Deltaca_theo(sec, self.tstop)
            (Deltav[mode],       Deltaca[mode])       = self.get_Deltav_Deltaca(sec, i)

            volbyarea[mode] = sec.diam/4
            vtocai[mode] = self.cm/(1E-1*2*h.FARADAY*volbyarea[mode])

            (diffv[mode], diffca[mode]) = self.get_diffv_diffca(i)
            print("----------------------------------------------------------------------")
            i = i+1

        return(Deltav, Deltaca, Deltav_theo, Deltaca_theo, volbyarea, vtocai, diffv, diffca)

    def test_startup(self):
        mechanism = 'caPump1'
        self.kappa = KappaNEURON.Kappa(membrane_species=[self.ca], kappa_file=os.path.dirname(KappaNEURON.__file__) + "/tests/" + mechanism + ".ka", regions=self.r)
        self.assertIsInstance(self.kappa._kappa_fluxes[0], KappaNEURON.KappaFlux)
        init()
        self.assertEqual(rxd.rxd._curr_indices, [1])

    def test_twoMembraneSpecies(self):
        self.napulse = h.NaPulse(self.sk(0.5))
        ## import pdb; pdb.set_trace()
        self.napulse.t0 = 0.5
        self.napulse.t1 = 1.2
        self.sk(0.5).t0_capulse = 1.0
        self.sk(0.5).t1_capulse = 1.5
        self.sk(0.5).fghk_capulse = 1
        self.sk(0.5).gbar_capulse = 0.0001
        self.na   = rxd.Species(self.r, name='na', charge=1, initial=0)
        self.kappa = KappaNEURON.Kappa(membrane_species=[self.ca, self.na], kappa_file=os.path.dirname(KappaNEURON.__file__) + "/tests/" + "napulse.ka", regions=self.r)
        self.assertIsInstance(self.kappa._kappa_fluxes[0], KappaNEURON.KappaFlux)
        self.assertIsInstance(self.kappa._kappa_fluxes[1], KappaNEURON.KappaFlux)
        
        ## Set up recordings
        self.rec_t = h.Vector()
        self.rec_t.record(h._ref_t)
        self.rec_cai = []
        self.rec_nai = []
        self.rec_v = []
        for sec in [self.sk]:
            self.rec_cai.append(h.Vector())
            self.rec_cai[-1].record(sec(0.5)._ref_cai)
            self.rec_nai.append(h.Vector())
            self.rec_nai[-1].record(sec(0.5)._ref_nai)
            self.rec_v.append(h.Vector())
            self.rec_v[-1].record(sec(0.5)._ref_v)

        init()
        ## Mapping from _curr_ptr_storage and _rxd_induced_currents
        ## onto state vector b in _rxd_reaction
        ## Depends on order of running tests
        ## self.assertEqual(rxd.rxd._curr_indices, [1, 7])
        ## Check length of get_memb_flux
        self.assertEqual(len(self.kappa._kappa_fluxes[0]._get_memb_flux(rxd.rxd._node_get_states())), 1)
        self.assertEqual(len(self.kappa._kappa_fluxes[1]._get_memb_flux(rxd.rxd._node_get_states())), 1)
        ## Check length of _cur_ptrs
        self.assertEqual(len(self.kappa._kappa_fluxes[0]._cur_ptrs), 1)
        self.assertEqual(len(self.kappa._kappa_fluxes[1]._cur_ptrs), 1)
        ## Check length of _cur_mapped
        self.assertEqual(len(self.kappa._kappa_fluxes[0]._cur_mapped), 1)
        self.assertEqual(len(self.kappa._kappa_fluxes[1]._cur_mapped), 1)
        ## Check length of first element of _cur_mapped is same as length of _cur_ptrs
        self.assertEqual(len(self.kappa._kappa_fluxes[0]._cur_mapped[0]), 
                         len(self.kappa._kappa_fluxes[0]._cur_ptrs[0]))
        self.assertEqual(len(self.kappa._kappa_fluxes[1]._cur_mapped[0]), 
                         len(self.kappa._kappa_fluxes[1]._cur_ptrs[0]))
        ## Check net charges are correct for Ca (2) and Na (1)
        self.assertEqual(self.kappa._kappa_fluxes[0]._net_charges, 2)
        self.assertEqual(self.kappa._kappa_fluxes[1]._net_charges, 1)

        run(3)
        times = np.array(self.rec_t)
        cai = numpy.array(self.rec_cai[0])
        nai = numpy.array(self.rec_nai[0])

        plt.subplots_adjust(left=0.25)
        nrow = 3
        fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(2.25*2, 2.5*2))
        ax[0].plot(self.rec_t, self.rec_v[0])
        ax[1].plot(self.rec_t, cai)
        ax[2].plot(self.rec_t, nai)
        fig.show()
        inds = np.where((times >= self.napulse.t0) & (times <= self.napulse.t1))
        self.assertGreater(max(nai[inds]), 0)
        inds = np.where((times >= self.sk.t0_capulse) & (times <= self.sk.t1_capulse))
        self.assertGreater(max(cai[inds]), 0)

    def test_twoMembraneSpeciesOneUncharged(self):
        """In this test we would like to give the Glutamate zero charge, so
        that there is no effect on the membrane potential.

        """ 
        self.glupulse = h.GluPulse(self.sk(0.5))
        ## import pdb; pdb.set_trace()
        self.glupulse.t0 = 0.5
        self.glupulse.t1 = 1.2
        self.sk(0.5).t0_capulse = 1.0
        self.sk(0.5).t1_capulse = 1.5
        self.sk(0.5).fghk_capulse = 1
        self.sk(0.5).gbar_capulse = 0.0001
        self.glu   = KappaNEURON.UnchargedSpecies(self.r, name='glu', initial=0)
        self.kappa = KappaNEURON.Kappa(membrane_species=[self.ca, self.glu], kappa_file=os.path.dirname(KappaNEURON.__file__) + "/tests/" + "nmda.ka", regions=self.r)
        self.assertIsInstance(self.kappa._kappa_fluxes[0], KappaNEURON.KappaFlux)
        # self.assertExists(self.kappa._kappa_fluxes[1])
        
        ## Set up recordings
        self.rec_t = h.Vector()
        self.rec_t.record(h._ref_t)
        self.rec_cai = []
        self.rec_glui = []
        self.rec_v = []
        for sec in [self.sk]:
            self.rec_cai.append(h.Vector())
            self.rec_cai[-1].record(sec(0.5)._ref_cai)
            self.rec_glui.append(h.Vector())
            self.rec_glui[-1].record(sec(0.5)._ref_glui)
            self.rec_v.append(h.Vector())
            self.rec_v[-1].record(sec(0.5)._ref_v)

        init()
        ## Mapping from _curr_ptr_storage and _rxd_induced_currents
        ## onto state vector b in _rxd_reaction
        ## Depends on order of running tests
        ## self.assertEqual(rxd.rxd._curr_indices, [1, 4])
        ## Check length of get_memb_flux
        self.assertEqual(self.kappa._kappa_fluxes[0]._membrane_flux, True)
        self.assertEqual(self.kappa._kappa_fluxes[1]._membrane_flux, False)
        ## Check membrane species are correct
        self.assertEqual(self.kappa._kappa_fluxes[0]._membrane_species[0], self.ca)
        self.assertEqual(self.kappa._kappa_fluxes[1]._membrane_species[0], self.glu)
        ## Check length of get_memb_flux
        self.assertEqual(len(self.kappa._kappa_fluxes[0]._get_memb_flux(rxd.rxd._node_get_states())), 1)
        ## Check length of _cur_ptrs
        self.assertEqual(len(self.kappa._kappa_fluxes[0]._cur_ptrs), 1)
        ## Check length of _cur_mapped
        self.assertEqual(len(self.kappa._kappa_fluxes[0]._cur_mapped), 1)
        ## Check length of first element of _cur_mapped is same as length of _cur_ptrs
        self.assertEqual(len(self.kappa._kappa_fluxes[0]._cur_mapped[0]), 
                         len(self.kappa._kappa_fluxes[0]._cur_ptrs[0]))
        ## Check net charges are correct for Ca (2) and Glu (1)
        self.assertEqual(self.kappa._kappa_fluxes[0]._net_charges, 2)
        self.kappa._kappa_fluxes[1]._net_charges = 1
        # self.assertEqual(self.kappa._kappa_fluxes[1]._net_charges, 0)

        run(3)
        # self.assertEqual(self.kappa._kappa_fluxes[1]._net_charges, 0)
        plt.subplots_adjust(left=0.25)
        nrow = 3
        fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(2.25*2, 2.5*2))
        ax[0].plot(self.rec_t, self.rec_v[0])
        ax[1].plot(self.rec_t, self.rec_cai[0])
        ax[2].plot(self.rec_t, self.rec_glui[0])
        ax[2].set_xlabel("time (ms)")
        ax[0].set_ylabel("V (mV)")
        ax[1].set_ylabel("cai")
        ax[2].set_ylabel("glui")

        fig.show()


    def test_injectCalcium(self):
        self.tstop = self.t1 + h.dt
        self.injectCalcium(ghk=0)
        self.do_plot()
        ## Needed for tests involving theoretical voltage and calcium to work
        self.tstop = self.t1
        ## Get data from both sections
        Deltav, Deltaca, Deltav_theo, Deltaca_theo, volbyarea, vtocai, diffv, diffca = self.get_stats()
        ## Check theory and simulation match, if using
        ## deterministic ('mod') simulation
        self.assertAlmostEqual(Deltav['mod'], Deltav_theo['mod'], 0)
        self.assertAlmostEqual(Deltaca['mod'], Deltaca_theo['mod'], 2)

        ## Calcium and voltage should be in sync
        self.assertAlmostEqual(Deltav['mod'], Deltaca['mod']/vtocai['mod'], 0)
        self.assertAlmostEqual(Deltav['kappa'], Deltaca['kappa']/vtocai['kappa'], 0)

        ## All calcium ion increments should be integers
        self.assertAlmostEqual(max(self.caitonum*diffca['kappa'] - np.round(self.caitonum*diffca['kappa'])), 0, 2)
        ## Calcium ion increments should be equal to voltage increments
        self.assertAlmostEqual(max(abs(vtocai['kappa']*diffv['kappa'][1:len(diffv['kappa'])-1] - diffca['kappa'][0:len(diffv['kappa'])-2])), 0, 2)


    def test_injectCalciumGHK(self):
        self.injectCalcium(ghk=1)
        self.do_plot()
        ## Get data from both sections
        Deltav, Deltaca, Deltav_theo, Deltaca_theo, volbyarea, vtocai, diffv, diffca = self.get_stats()

        ## Calcium and voltage should be in sync
        for mode in ['mod', 'kappa']:
            self.assertAlmostEqual(Deltav[mode], Deltaca[mode]/vtocai[mode], 0)

        ## All calcium ion increments should be integers
        self.assertAlmostEqual(max(self.caitonum*diffca['kappa'] - np.round(self.caitonum*diffca['kappa'])), 0, 2)
        ## Calcium ion increments should be equal to voltage increments
        self.assertAlmostEqual(max(abs(vtocai['kappa']*diffv['kappa'][1:len(diffv['kappa'])-1] - diffca['kappa'][0:len(diffv['kappa'])-2])), 0, 2)

    def test_injectCalciumPump(self):
        self.t1 = 2
        self.tstop = 2
        self.k1 = 1
        self.injectCalcium(ghk=0)
        self.do_plot()
        ## Get data from both sections
        Deltav, Deltaca, Deltav_theo, Deltaca_theo, volbyarea, vtocai, diffv, diffca = self.get_stats()

        ## Check theory and simulation match, if using
        ## deterministic ('mod') simulation
        self.assertAlmostEqual(Deltav['mod'],  Deltav_theo['mod'],  0)
        self.assertAlmostEqual(Deltaca['mod'], Deltaca_theo['mod'], 2)

        ## Calcium and voltage should be in sync
        for mode in ['mod', 'kappa']:
            self.assertEqualWithinTol(Deltav[mode], Deltaca[mode]/vtocai[mode])

        ## Check that kappa and deterministic simulations agree to
        ## within 10%
        tol = 0.15
        self.assertLess(abs((Deltav['kappa'] - Deltav['mod'])/(Deltav['mod'] - self.v0)), tol)
        self.assertLess(abs((Deltaca['kappa'] - Deltaca['mod'])/Deltaca['mod']), tol)
        

    def test_injectCalciumPumpGHK(self):
        self.t1 = 2
        self.tstop = 2
        self.k1 = 1
        self.injectCalcium(ghk=1)
        self.do_plot()
        ## Get data from both sections
        Deltav, Deltaca, Deltav_theo, Deltaca_theo, volbyarea, vtocai, diffv, diffca = self.get_stats()

        ## Calcium and voltage should be in sync
        for mode in ['mod', 'kappa']:
            self.assertAlmostEqual(Deltav[mode], Deltaca[mode]/vtocai[mode], 0)

        ## Check that kappa and deterministic simulations agree to
        ## within 10%
        tol = 0.1
        self.assertLess(abs((Deltav['kappa'] - Deltav['mod'])/(Deltav['mod'] - self.v0)), 0.1)
        self.assertLess(abs((Deltaca['kappa'] - Deltaca['mod'])/Deltaca['mod']), 0.1)

        ## All calcium ion increments should be integers
        self.assertAlmostEqual(max(self.caitonum*diffca['kappa'] - np.round(self.caitonum*diffca['kappa'])), 0, 2)
        ## Calcium ion increments should be equal to voltage increments
        self.assertAlmostEqual(max(abs(vtocai['kappa']*diffv['kappa'][1:len(diffv['kappa'])-1] - diffca['kappa'][0:len(diffv['kappa'])-2])), 0, 2)


    def test_injectCalciumPump2(self):
        self.t1 = 2.0
        self.tstop = 3.0
        self.k1 = 1
        self.k2 = 0
        self.P0 = 0.20 
        self.injectCalcium(ghk=0, mechanism='caPump2')
        self.do_plot()

        ## Run through both sections
        times = np.array(self.rec_t)
        i = 0
        for sec in h.allsec():
            mode = self.get_mode(sec)
            print(mode)
            v = np.array(self.rec_v[i])
            self.assertAlmostEqual(v[np.where(np.isclose(times, self.t1 + 0.1))],
                                   v[np.where(np.isclose(times, self.tstop))])
            cai = np.array(self.rec_cai[i])
            self.assertGreater(cai[np.where(np.isclose(times, self.t1 + 0.1))],
                               cai[np.where(np.isclose(times, self.tstop))])
            i += 1

    def test_injectCalciumPump2k2(self):
        self.t1 = 2.0
        self.tstop = 5.0
        self.k1 = 1
        self.k2 = 10
        self.P0 = 0.20 
        self.injectCalcium(ghk=0, mechanism='caPump2')
        self.do_plot()
        ## Run through both sections
        times = np.array(self.rec_t)
        i = 0
        for sec in h.allsec():
            mode = self.get_mode(sec)
            print(mode)
            v = np.array(self.rec_v[i])
            self.assertGreater(v[np.where(np.isclose(times, self.t1 + 0.1))],
                                   v[np.where(np.isclose(times, self.tstop))])
            cai = np.array(self.rec_cai[i])
            self.assertGreater(cai[np.where(np.isclose(times, self.t1 + 0.1))],
                               cai[np.where(np.isclose(times, self.tstop))])
            i += 1

    def tearDown(self):
        self.kappa.__del__()
        print(len(KappaNEURON._kappa_schemes))
        if not self.mechanism is None:
            h('uninsert ' + self.mechanism, sec=self.sm)
        self.mechanism = None


testSuite = unittest.TestSuite()
testSuite.addTest(TestCaAccumulation('test_injectCalcium'))
# testSuite.addTest(TestCaAccumulation('test_injectCalciumGHK'))
# testSuite.addTest(TestCaAccumulation('test_injectCalciumPump'))
# testSuite.addTest(TestCaAccumulation('test_injectCalciumPumpGHK'))
# testSuite.addTest(TestCaAccumulation('test_injectCalciumPump2'))
# testSuite.addTest(TestCaAccumulation('test_injectCalciumPump2k2'))
# testSuite.addTest(TestCaAccumulation('test_twoMembraneSpecies'))
# testSuite.addTest(TestCaAccumulation('test_twoMembraneSpeciesOneUncharged'))

unittest.TextTestRunner(verbosity=2).run(testSuite)

# def quick():
#     print ("hello")
#     unittest.TextTestRunner(verbosity=2).run(testSuite)

# if __name__ == '__main__':
#     unittest.main()
