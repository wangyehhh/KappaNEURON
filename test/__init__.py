import KappaNEURON
import os
import unittest
import neuron
from neuron import *
from neuron import rxd
import numpy as np
import matplotlib.pyplot as plt

class TestCaAccumulation(unittest.TestCase):
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
    ca = rxd.Species(r, name='ca', charge=2, initial=0.0)

    ## Set up a mod simulation of one compartment
    sm = h.Section()
    sm.insert("capulse")
    sm.L=0.2
    sm.diam=0.5
    ## Calcium accumulation via a pump with the pumping turned off
    sm.insert("caPump1")
    
    ## Time of pulse
    t0 = 1.0
    t1 = 1.1
    
    ## Conversions
    NA =  6.02214129e23     # Avogadro's number

    gbar = 0.001
    cm = sk(0.5).cm

    def setUp(self):
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
        self.caitonum = self.NA*np.pi*(self.sk.diam**2)/4*self.sk.L*1e-18 

    def injectCalcium(self, ghk=0, k1=0):
        self.kappa = KappaNEURON.Kappa([self.ca], self.__module__ + "/caPump1.ka", self.r, verbose=True)
        self.kappa.setVariable('k1', k1)
        self.sm(0.5).k1_caPump1 = k1
        
        KappaNEURON.progress = False
        KappaNEURON.verbose = False

        self.assertIsInstance(self.sk, nrn.Section)
        self.assertEqual(self.ca.initial, 0.0)
        for sec in h.allsec():
            ## This forces eca to be a constant, rather than being
            ## computed from Nernst equation at every time step
            ## See http://www.neuron.yale.edu/neuron/static/new_doc/modelspec/programmatic/ions.html?highlight=ion_style#ion_style
            if (ghk == 1):
                sec.fghk_capulse = 1
            else:
                sec.push(); h('ion_style("ca_ion", 3, 1, 0, 0, 1)') ; h.pop_section()
                
            for seg in sec:
                seg.t0_capulse = self.t0
                seg.t1_capulse = self.t1

                seg.gbar_capulse = self.gbar
        h.dt = h.dt/20

        ## Initialise simulation
        init()
        self.v0 = self.sk(0.5).v
        self.assertEqual(h.t, 0.0)
        self.assertEqual(self.sk(0.5).cai, 0.0)
        ## Run
        run(1.15)
        self.assertAlmostEqual(h.t, 1.15)
        self.assertGreater(self.sk(0.5).cai, 0.0)


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
        fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(2.25*2, 2.5*2))
        i = 0
        for sec in h.allsec():
            (diffv, diffca) = self.get_diffv_diffca(i)
            ax[0].plot(self.rec_t, self.rec_v[i], color='br'[i])
            ax[1].plot(self.rec_t, self.rec_cai[i], color='br'[i])
            ax[2].plot(diffv[1:len(diffv)-1], self.caitonum*diffca[0:len(diffv)-2], 'o', color='br'[i])
            fig.show()        
            i = i + 1

    def get_Deltav_Deltaca_theo(self, sec):
        eca = sec(0.5).eca
        volbyarea = sec.diam/4
        vtocai = self.cm/(1E-1*2*h.FARADAY*volbyarea)
        print 'Theoretical voltage and Ca difference:'
        Deltav_theo = (eca - self.v0)*(1 - np.exp(-(self.t1 - self.t0)*1000*self.gbar/self.cm))
        Deltaca_theo = Deltav_theo*vtocai
        print Deltav_theo, Deltaca_theo
        return(Deltav_theo, Deltaca_theo)

    def get_Deltav_Deltaca(self, sec, i):
        v1 = sec(0.5).v
        print 'Actual voltage and Ca difference:'
        Deltav = v1 - self.v0
        Deltaca = sec(0.5).cai - self.rec_cai[i][0]
        print Deltav, Deltaca
        return(Deltav, Deltaca)

    def test_injectCalcium(self):
        self.injectCalcium(ghk=0)

        if self.plot:
            self.do_plot()

        ## Run through both sections
        i = 0
        for sec in h.allsec():
            ## Determine if section contains mod pump or kappa pump
            mode = 'kappa'
            for mech in sec(0.5):
                if mech.name() == 'caPump1':
                    mode = 'mod'
            print mode

            ## Print some variables
            v1 = sec(0.5).v
            print("Eca=%f, t0=%f, t1=%f, gbar=%f, cm=%f, v0=%f, v1=%f" % (sec(0.5).eca, self.t0, self.t1, self.gbar, self.cm, self.v0, v1))

            ## Check theory and simulation match, if using
            ## deterministic ('mod') simulation
            (Deltav_theo, Deltaca_theo) = self.get_Deltav_Deltaca_theo(sec)
            (Deltav,      Deltaca     ) = self.get_Deltav_Deltaca(sec, i)
            if mode == 'mod':
                self.assertAlmostEqual(Deltav, Deltav_theo, 0)
                self.assertAlmostEqual(Deltaca, Deltaca_theo, 2)

            ## Check voltage and calcium are in sync
            volbyarea = sec.diam/4
            vtocai = self.cm/(1E-1*2*h.FARADAY*volbyarea)
            self.assertAlmostEqual(Deltav, Deltaca/vtocai, 1)

            ## Check all differences are the same
            (diffv, diffca) = self.get_diffv_diffca(i)
            if mode == 'kappa':
                ## All calcium ion increments should be integers
                self.assertAlmostEqual(max(self.caitonum*diffca - np.round(self.caitonum*diffca)), 0, 2)
                ## Calcium ion increments should be equal to voltage increments
                self.assertAlmostEqual(max(abs(vtocai*diffv[1:len(diffv)-1] - diffca[0:len(diffv)-2])), 0, 2)
            i = i + 1

    def test_injectCalciumGHK(self):
        self.injectCalcium(ghk=1)
        ## Run through both sections
        i = 0
        for sec in h.allsec():
            ## Determine if section contains mod pump or kappa pump
            mode = 'kappa'
            for mech in sec(0.5):
                if mech.name() == 'caPump1':
                    mode = 'mod'
            print mode
            
            ## Print some variables
            v1 = sec(0.5).v
            print("Eca=%f, t0=%f, t1=%f, gbar=%f, cm=%f, v0=%f, v1=%f" % (sec(0.5).eca, self.t0, self.t1, self.gbar, self.cm, self.v0, v1))

            ## The theoretical values here don't apply since we are
            ## using GHK, but get them anyway for information
            (Deltav_theo, Deltaca_theo) = self.get_Deltav_Deltaca_theo(sec)
            (Deltav,      Deltaca     ) = self.get_Deltav_Deltaca(sec, i)

            ## Check voltage and calcium are in sync
            volbyarea = sec.diam/4
            vtocai = self.cm/(1E-1*2*h.FARADAY*volbyarea)
            self.assertAlmostEqual(Deltav, Deltaca/vtocai, 1)

            ## Check all differences are the same
            (diffv, diffca) = self.get_diffv_diffca(i)
            if mode == 'kappa':
                ## All calcium ion increments should be integers
                self.assertAlmostEqual(max(self.caitonum*diffca - np.round(self.caitonum*diffca)), 0, 2)
                ## Calcium ion increments should be equal to voltage increments
                self.assertAlmostEqual(max(abs(vtocai*diffv[1:len(diffv)-1] - diffca[0:len(diffv)-2])), 0, 2)
            i = i + 1

    def test_injectCalciumPump(self):
        self.injectCalcium(ghk=0, k1=0.01)


    @unittest.skip("skip for now")
    def test_injectCalcium2(self):
        init()
        self.assertEqual(h.t, 0.0)
        run(3)
        self.assertAlmostEqual(h.t, 3.0)

    def tearDown(self):
        self.kappa = None

testSuite = unittest.TestSuite()
#testSuite.addTest(TestCaAccumulation('test_injectCalcium'))
testSuite.addTest(TestCaAccumulation('test_injectCalciumGHK'))
        
if __name__ == '__main__':
    unittest.main()
