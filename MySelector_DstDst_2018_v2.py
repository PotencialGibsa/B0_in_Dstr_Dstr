from ROOT import *
#from variables import *
import glob
import numpy as np
import math

def sqrt(var):
      return math.sqrt(abs(var))
#############################################################
PDG_MUON_MASS    =   0.1056583745
PDG_PION_MASS    =   0.13957061
PDG_PIOZ_MASS    =   0.1349770
PDG_KAON_MASS    =   0.493677
PDG_PROTON_MASS  =   0.9382720813
PDG_KSHORT_MASS  =   0.497611
PDG_KSHORT_DM    =   0.000013
PDG_KSHORT_TIME  =   0.8954 * 0.0000000001
PDG_KS_MASS      =   PDG_KSHORT_MASS
PDG_LAMBDA_MASS  =   1.115683
PDG_LAMBDA_DM    =   0.000006
PDG_LAMBDA_TIME  =   2.632 * 0.0000000001
PDG_SIGMA0_MASS  =   1.192642
PDG_XImunus_MASS =   1.32171
PDG_XImunus_DM   =   0.00007
PDG_XImunus_TIME =   1.639 * 0.0000000001
PDG_OMmunus_MASS =   1.67245
PDG_OMmunus_DM   =   0.00029
PDG_OMmunus_TIME =   0.821 * 0.0000000001
PDG_DPM_MASS     =   1.86965
PDG_DPM_DM       =   0.00005
PDG_DPM_TIME     =   1.040 * 0.000000000001
PDG_DZ_MASS      =   1.86483
PDG_DZ_DM        =   0.00005
PDG_DZ_TIME      =   0.4101 * 0.000000000001
PDG_DS_MASS      =   1.96834
PDG_DS_DM        =   0.00007
PDG_DS_TIME      =   0.504 * 0.000000000001
PDG_LAMCP_MASS   =   2.28646
PDG_LAMCP_DM     =   0.00031
PDG_LAMCP_TIME   =   2.00 * 0.0000000000001
PDG_XICZ_MASS    =   2.47087
PDG_XICZ_DM      =   0.00031
PDG_XICZ_TIME    =   1.12 * 0.0000000000001
PDG_XICP_MASS    =   2.46787
PDG_XICP_DM      =   0.00030
PDG_XICP_TIME    =   4.42 * 0.0000000000001
PDG_KSTARZ_MASS  =   0.89555
PDG_KSTARZ_GAMMA =   0.0473
PDG_KSTARP_MASS  =   0.89176
PDG_KSTARP_GAMMA =   0.0503
PDG_PHI_MASS     =   1.019461
PDG_PHI_GAMMA    =   0.004249
PDG_JPSI_MASS    =   3.096900
PDG_PSI2S_MASS   =   3.686097
PDG_X3872_MASS   =   3.87169
PDG_BU_MASS      =   5.27932
PDG_BU_TIME      =   1.638 * 0.000000000001
PDG_B0_MASS      =   5.27963
PDG_B0_TIME      =   1.520 * 0.000000000001
PDG_BS_MASS      =   5.36689
PDG_BS_TIME      =   1.509 * 0.000000000001
PDG_BC_MASS      =   6.2749
PDG_BC_TIME      =   0.507 * 0.000000000001
PDG_LB_MASS      =   5.61960
PDG_LB_TIME      =   1.470 * 0.000000000001
PDG_XIBZ_MASS    =   5.7919
PDG_XIBZ_TIME    =   1.479 * 0.000000000001
PDG_XIBM_MASS    =   5.7970
PDG_XIBM_TIME    =   1.571 * 0.000000000001
PDG_OMBM_MASS    =   6.0461
PDG_OMBM_TIME    =   1.64 * 0.000000000001
PDG_C            =   29979245800. ### in cm/c
PDG_DSTR         =   2.01026
###  }}}

md      = RooRealVar ( "md"     ,"M(D) [GeV]"               , PDG_DZ_MASS-0.03 , PDG_DZ_MASS+0.03  )
mds     = RooRealVar ( "mds"    ,"M(Dstar) [GeV]"           , 2.004, 2.019 )
#
dspt    = RooRealVar ( "dspt"   ,"dspt"                     , 3.0   , 33.0 )
dzpt    = RooRealVar ( "dzpt"   ,"dzpt"                     , 2.5   , 22.5  )
kspt    = RooRealVar ( "kspt"   ,"kspt"                     , 0.0   , 10.0 )
pipt    = RooRealVar ( "pipt"   ,"pipt"                     , 0.1   , 1.5 )
#
dset    = RooRealVar ( "dset"   ,"dset"                     , -2.6  , 2.6   )
dzet    = RooRealVar ( "dzet"   ,"dzet"                     , -2.6  , 2.6   )
kset    = RooRealVar ( "kset"   ,"kset"                     , -2.6  , 2.6   )
piet    = RooRealVar ( "piet"   ,"piet"                     , -2.6  , 2.6   )
#
dzds2   = RooRealVar ( "dzds2"  ,"DetSign"                  , 0.0   , 50    )
dzds3   = RooRealVar ( "dzds3"  ,"DetSign"                  , 0.0   , 80    )
#
dsvtxp  = RooRealVar ( "dsvtxp" ,"dsvtxp"                   , -0.1  , 1.1   )
dzvtxp  = RooRealVar ( "dzvtxp" ,"dzvtxp"                   , -0.1  , 1.1   )
#
ipsmin  = RooRealVar ( "ipsmin" ,"ipsmin"                   , 0.0   , 1000.0)
vari    = RooRealVar ( "vari"   ,"vari"                     , 0.0   , 5)
vdr     = RooRealVar ( "vdr"    ,"vdr"                      , 0.0   , 2)
vdz     = RooRealVar ( "vdz"    ,"vdz"                      , 0.0   , 5)

def DetachSignificance2(vtx, vtxE1, vtxE2):
    return sqrt( vtx.X()**2 / (vtxE1.X()**2 + vtxE2.X()**2) + vtx.Y()**2 / (vtxE1.Y()**2 + vtxE2.Y()**2))

def DetachSignificance3(vtx, vtxE1, vtxE2):
    return sqrt( vtx.X()**2 / (vtxE1.X()**2 + vtxE2.X()**2) + vtx.Y()**2 / (vtxE1.Y()**2 + vtxE2.Y()**2) + vtx.Z()**2 / (vtxE1.Z()**2 + vtxE2.Z()**2))

def DirectionCos2 (v1, v2):
    r1 = sqrt(v1.X()**2 + v1.Y()**2)
    r2 = sqrt(v2.X()**2 + v2.Y()**2)
    return ( v1.X() * v2.X() + v1.Y() * v2.Y() ) / (r1*r2 + 0.0000001)

def DirectionCos3 (v1, v2):
    r1 = sqrt(v1.X()**2 + v1.Y()**2 + v1.Z()**2)
    r2 = sqrt(v2.X()**2 + v2.Y()**2 + v2.Z()**2)
    return ( v1.X() * v2.X() + v1.Y() * v2.Y() + v1.Z() * v2.Z() ) / (r1*r2 + 0.0000001)

def DirectionChi22 (vtx0, vtx0E, vtx1, vtx1E, P, PE):
    dvtx    = vtx1 - vtx0
    dvtxE   = TVector3( sqrt(vtx0E.X()**2 + vtx1E.X()**2), sqrt(vtx0E.Y()**2 + vtx1E.Y()**2), sqrt(vtx0E.Z()**2 + vtx1E.Z()**2) )
    Pscaled = P * (dvtx.Mag() / P.Mag())
    PscaledE= PE * (dvtx.Mag() / P.Mag())
    return DetachSignificance2 (Pscaled - dvtx, PscaledE, dvtxE)

def DirectionChi23 (vtx0, vtx0E, vtx1, vtx1E, P, PE):
    dvtx    = vtx1 - vtx0 ## vertex difference
    dvtxE   = TVector3( sqrt(vtx0E.X()**2 + vtx1E.X()**2), sqrt(vtx0E.Y()**2 + vtx1E.Y()**2), sqrt(vtx0E.Z()**2 + vtx1E.Z()**2) ) ## its error
    Pscaled = P * (dvtx.Mag() / P.Mag()) ## scaled momentum to be the same length as vertex difference
    PscaledE= PE * (dvtx.Mag() / P.Mag()) ## its error
    return DetachSignificance3 (Pscaled - dvtx, PscaledE, dvtxE)


def TH1NormBinWidth(hh):
    _htemp = hh.Clone()
    _int = _htemp.Integral()
    for _i in range(1, hh.GetNbinsX()+1):
        _htemp[_i] = _htemp[_i] / _htemp.GetBinWidth(_i)
    #
    _nam = hh.GetName() + '_no'
    _htemp.SetName('_nam')
    _htemp.Scale(_int / _htemp.Integral())
    return _htemp

def stri(n):
    if n>=0 and n<10:
        return '0'+str(n)
    else:
        return str(n)

Parts = ['A','B','C','D1','D2','D3','D4','D5'] #for ParkingBPH

for iter in Parts:
    Parts1 = glob.glob('/eos/home-d/dshmygol/crab_projects/2018_DstDstv1_%s/*/*/*/000?' % iter )

    if 'D' in iter:
        print(iter)
        Parts1 = glob.glob('/eos/home-d/dshmygol/crab_projects/2018_DstDstv1_D/ParkingBPH%s/*/*/000?' % iter[1] )

    #without it doesnt work

    for iter1 in Parts1:
        MyFileNames = glob.glob(iter1 + '/BFinder*.root')
        print('The number of files is %i for %s' % (len(MyFileNames), iter))
        print('\nAdding the files to the chain')
        chain = TChain('wztree')
        for name in MyFileNames:
            chain.Add(name)
        #
        strPark = 'DstDstv2'+ '_' + iter
        _fileOUT = '/eos/home-d/dshmygol/DstDstv2/SimpleFile_%s.root' % (strPark)
        fileOUT = TFile (_fileOUT, "recreate")
        mytree = TTree("mytree","mytree")

        nevt = chain.GetEntries()
        print('Number of events: %s' % (nevt))

        if 1>0:
            # declaring variables {{{
            my_vars = [
                    "SAMEEVENT",
                    #Mass D*
                    'Dstr1_masd', 'Dstr1_masd_CV',
                    'Dstr2_masd', 'Dstr2_masd_CV',

                    'KA1_pt',  'KA1_ips', # 'KA1_eta',
                    'PI1_pt', 'PI1_ips',  'PI1_dz',
                    'PI2_pt', 'PI2_ips',  'PI2_dz',
                    'KA2_pt', 'KA2_ips', 'KA2_chrg', # 'KA2_eta',
                    'PIS1_pt', 'PIS1_eta',
                    'PIS1_ips', 'PIS1_dr', 'PIS1_dz', 'PIS1_chrg',
                    'PIS2_pt', 'PIS2_eta',
                    'PIS2_ips', 'PIS2_dr', 'PIS2_dz', 'PIS2_chrg',

                    'D_mass0', 'D_mass', 'D_mass_CV',
                    'D_pt',  'D_eta', # 'D_vtxErrXY', 'D_vtxErrZ',
                    'D_pvdistsignif2', 'D_pvdistsignif3', 'D_BSdistsignif2',
                    'D_pvdist', 'D_vtxprob',
                    'D_pvcos2', 'D_pvcos3', 'D_BScos2',

                    'E_mass0', 'E_mass', 'E_mass_CV',
                    'E_pt', 'E_eta', # 'E_vtxErrXY', 'E_vtxErrZ',
                    'E_pvdistsignif2', 'E_pvdistsignif3', 'E_BSdistsignif2',
                    'E_pvdist', 'E_vtxprob',
                    'E_pvcos2', 'E_pvcos3', 'E_BScos2',
            #         'Dstr_mass', 'Dstr_mass_delta', 'Dstr_mass_delta1',

                    'B_mass',  'B_masd', 'B_masdd',
                    'B_pt', 'B_eta', 'B_p',
                    'B_pvdistsignif2', 'B_pvdistsignif3',
                    'B_pvdist', 'B_vtxprob',
                    'B_pvcos2', 'B_pvcos3',

                    'D_Bdistsignif2', 'D_Bdistsignif3',
                    'D_Bcos2', 'D_Bcos3',

                    'E_Bdistsignif2', 'E_Bdistsignif3',
                    'E_Bcos2', 'E_Bcos3',

                    'TRIG_0', 'TRIG_1', 'TRIG_2', 'TRIG_3', 'TRIG_4', 'TRIG_5', 'TRIG_6', 'TRIG_7', 'TRIG_8', 'TRIG_9',
                    'GEN_minDRsum',
            ]

            for _var in my_vars:
                exec(_var+'=np.zeros(1, dtype=float)')
                exec('mytree.Branch("' + _var + '"' + ' '*(25-len(_var)) + ',' + _var + ' '*(25-len(_var)) + ', "'+ _var + '/D")')

            ### finished with deploying vars for SimpleFile }}}

            # P4 ROOT Lorentz vector
            KA1P4, KA1P4_CV, KA2P4, KA2P4_CV, PI1P4, PI1P4_CV, PI2P4, PI2P4_CV, PIS1P4, PIS1P4_CV,PIS2P4, PIS2P4_CV, DP4, DP4_CV,EP4, EP4_CV, BP4 = [TLorentzVector() for i in range(17)]
            # Trigger fire counters
            nEvt_trig0, nEvt_trig1, nEvt_trig2, nEvt_trig3, nEvt_trig4, nEvt_trig5, nEvt_trig6, nEvt_trig7, nEvt_trig8, nEvt_trig9 = [int(0) for i in range(10)]
            nEvt_trig0_after, nEvt_trig1_after, nEvt_trig2_after, nEvt_trig3_after, nEvt_trig4_after, nEvt_trig5_after, nEvt_trig6_after, nEvt_trig7_after, nEvt_trig8_after, nEvt_trig9_after = [int(0) for i in range(10)]

            BBB = 0

            # looping all over the events
            for evt in range(nevt):
                # breaking if the event has no entries
                if chain.GetEntry(evt) <= 0:
                    break
                # securing that _nCand is set to 0 for every event before assign the real chain value
                _nCand = 0
                _nCand = chain.nCand
                # if no Dstr candidate is storaged, then move to next event
                if _nCand < 1:
                    continue
                for cand in range(_nCand):
                    #it only works when we set variable with [0] by its side
                    # assigning Triggers variables
                    TRIG_0[0] = 1 if chain.trig0_fi[cand] else 0
                    TRIG_1[0] = 1 if chain.trig1_fi[cand] else 0
                    TRIG_2[0] = 1 if chain.trig2_fi[cand] else 0
                    TRIG_3[0] = 1 if chain.trig3_fi[cand] else 0
                    TRIG_4[0] = 1 if chain.trig4_fi[cand] else 0
                    TRIG_5[0] = 1 if chain.trig5_fi[cand] else 0
                    TRIG_6[0] = 1 if chain.trig6_fi[cand] else 0
                    TRIG_7[0] = 1 if chain.trig7_fi[cand] else 0
                    TRIG_8[0] = 1 if chain.trig8_fi[cand] else 0
                    TRIG_9[0] = 1 if chain.trig9_fi[cand] else 0
                    # assigning PV info
                    PV = TVector3(chain.PV_becos_XX[cand], chain.PV_becos_YY[cand], chain.PV_becos_ZZ[cand])
                    PVE = TVector3(sqrt(chain.PV_becos_EX[cand]), sqrt(chain.PV_becos_EY[cand]), sqrt(chain.PV_becos_EZ[cand]))
                    BESP = TVector3(chain.BESP_x, chain.BESP_y, chain.BESP_z)
                    BESPE = TVector3(chain.BESP_ex, chain.BESP_ey, chain.BESP_ez)

                    # assigning info about the Pi and KA that comes from D0

                    # assigning KA1 info
                    KA1P4.SetXYZM (chain.KA1_px[cand], chain.KA1_py[cand], chain.KA1_pz[cand], PDG_PION_MASS)
                    KA1P4_CV.SetXYZM (chain.KA1_px_CV[cand], chain.KA1_py_CV[cand], chain.KA1_pz_CV[cand], PDG_PION_MASS)
                    KA1_pt[0] = KA1P4.Pt()
                    # KA1_eta[0] = KA1P4.Eta()
                    KA1_ips[0] = chain.KA1_ips[cand]

                    # assigning KA2 info
                    KA2P4.SetXYZM (chain.KA3_px[cand], chain.KA3_py[cand], chain.KA3_pz[cand], PDG_PION_MASS)
                    KA2P4_CV.SetXYZM (chain.KA3_px_CV[cand], chain.KA3_py_CV[cand], chain.KA3_pz_CV[cand], PDG_PION_MASS)
                    KA2_pt[0] = KA2P4.Pt()
                    # KA2_eta[0] = KA2P4.Eta()
                    KA2_ips[0] = chain.KA3_ips[cand]
                    KA2_chrg[0] = chain.KA3_charge[cand]

                    # assigning PI1 info
                    PI1P4.SetXYZM (chain.PI2_px[cand], chain.PI2_py[cand], chain.PI2_pz[cand], PDG_PION_MASS)
                    PI1P4_CV.SetXYZM (chain.PI2_px_CV[cand], chain.PI2_py_CV[cand], chain.PI2_pz_CV[cand], PDG_PION_MASS)
                    PI1_pt[0] = PI1P4.Pt()
                    PI1_ips[0] = chain.PI2_ips[cand]
                    PI1_dz[0] = chain.PI2_dz[cand]

                    # assigning PI2 info
                    PI2P4.SetXYZM (chain.PI4_px[cand], chain.PI4_py[cand], chain.PI4_pz[cand], PDG_PION_MASS)
                    PI2P4_CV.SetXYZM (chain.PI4_px_CV[cand], chain.PI4_py_CV[cand], chain.PI4_pz_CV[cand], PDG_PION_MASS)
                    PI2_pt[0] = PI2P4.Pt()
                    PI2_ips[0] = chain.PI4_ips[cand]
                    PI2_dz[0] = chain.PI4_dz[cand]


                    # assigning info about D0
                    DV = TVector3(chain.D_DecayVtxX[cand], chain.D_DecayVtxY[cand], chain.D_DecayVtxZ[cand])
                    DVE = TVector3(sqrt(chain.D_DecayVtxXE[cand]), sqrt(chain.D_DecayVtxYE[cand]), sqrt(chain.D_DecayVtxZE[cand]))
                    DP4.SetXYZM( chain.D_px[cand], chain.D_py[cand], chain.D_pz[cand], chain.D_mass[cand])
                    DP4_CV.SetXYZM( chain.D_px_CV[cand], chain.D_py_CV[cand], chain.D_pz_CV[cand], chain.D_mass_CV[cand])
                    DP3 = DP4_CV.Vect()
                    D_pvdistsignif2[0] = DetachSignificance2( DV - PV, PVE, DVE)
                    D_pvdistsignif3[0] = DetachSignificance3( DV - PV, PVE, DVE)
                    D_BSdistsignif2[0] = DetachSignificance2( DV - BESP, BESPE, DVE) ## wrt beamspot
                    D_pvcos2[0] = DirectionCos2 (DV - PV, DP3)
                    D_pvcos3[0] = DirectionCos3 (DV - PV, DP3)
                    D_BScos2[0] = DirectionCos2 (DV - BESP, DP3) ## wrt beamspot
                    D_pvdist[0] = (DV - PV).Mag()
                    D_pt[0] = DP4.Pt()
                    D_eta[0] = DP4.Eta()
                    D_vtxprob[0]= chain.D_Prob[cand]
                    D_mass[0] = chain.D_mass[cand]
                    D_mass0[0] = (KA1P4 + PI1P4).M()
                    D_mass_CV[0] = chain.D_mass_CV[cand]
                    if D_pvdistsignif2[0] < 3 : continue
                    if D_pvcos2[0] < 0.9 : continue

                    # assigning info about E0
                    EV = TVector3(chain.E_DecayVtxX[cand], chain.E_DecayVtxY[cand], chain.E_DecayVtxZ[cand])
                    EVE = TVector3(sqrt(chain.E_DecayVtxXE[cand]), sqrt(chain.E_DecayVtxYE[cand]), sqrt(chain.E_DecayVtxZE[cand]))
                    EP4.SetXYZM( chain.E_px[cand], chain.E_py[cand], chain.E_pz[cand], chain.E_mass[cand])
                    EP4_CV.SetXYZM( chain.E_px_CV[cand], chain.E_py_CV[cand], chain.E_pz_CV[cand], chain.E_mass_CV[cand])
                    EP3 = EP4_CV.Vect()
                    E_pvdistsignif2[0] = DetachSignificance2( EV - PV, PVE, EVE)
                    E_pvdistsignif3[0] = DetachSignificance3( EV - PV, PVE, EVE)
                    E_BSdistsignif2[0] = DetachSignificance2( EV - BESP, BESPE, EVE) ## wrt beamspot
                    E_pvcos2[0] = DirectionCos2 (EV - PV, EP3)
                    E_pvcos3[0] = DirectionCos3 (EV - PV, EP3)
                    E_BScos2[0] = DirectionCos2 (EV - BESP, EP3) ## wrt beamspot
                    E_pvdist[0] = (EV - PV).Mag()
                    E_pt[0] = EP4.Pt()
                    E_eta[0] = EP4.Eta()
                    E_vtxprob[0]= chain.E_Prob[cand]
                    E_mass[0] = chain.E_mass[cand]
                    E_mass0[0] = (KA2P4 + PI2P4).M()
                    E_mass_CV[0] = chain.E_mass_CV[cand]
                    if E_pvdistsignif2[0] < 3 : continue
                    if E_pvcos2[0] < 0.9 : continue

                    # assigning PIS1 info
                    PIS1P4.SetXYZM (chain.PIS1_px[cand], chain.PIS1_py[cand], chain.PIS1_pz[cand], PDG_PION_MASS)
                    PIS1P4_CV.SetXYZM (chain.PIS1_px_CV[cand], chain.PIS1_py_CV[cand], chain.PIS1_pz_CV[cand], PDG_PION_MASS)
                    PIS1_pt[0]   = PIS1P4.Pt()
                    PIS1_eta[0]  = PIS1P4.Eta()
                    PIS1_ips[0]  = chain.PIS1_ips[cand]
                    PIS1_dr[0]   = chain.PIS1_dr[cand] ## wrt KS vertex -- not much useful probably
                    PIS1_dz[0]   = chain.PIS1_dz[cand] ## wrt KS vertex -- not much useful probably
                    PIS1_chrg[0] = chain.PIS1_charg[cand]

                    # assigning PIS2 info
                    PIS2P4.SetXYZM (chain.PIS2_px[cand], chain.PIS2_py[cand], chain.PIS2_pz[cand], PDG_PION_MASS)
                    PIS2P4_CV.SetXYZM (chain.PIS2_px_CV[cand], chain.PIS2_py_CV[cand], chain.PIS2_pz_CV[cand], PDG_PION_MASS)
                    PIS2_pt[0]   = PIS2P4.Pt()
                    PIS2_eta[0]  = PIS2P4.Eta()
                    PIS2_ips[0]  = chain.PIS2_ips[cand]
                    PIS2_dr[0]   = chain.PIS2_dr[cand] ## wrt KS vertex -- not much useful probably
                    PIS2_dz[0]   = chain.PIS2_dz[cand] ## wrt KS vertex -- not much useful probably
                    PIS2_chrg[0] = chain.PIS2_charg[cand]

                    Dstr1_masd[0] = (DP4 + PIS1P4).M() - DP4.M() + PDG_DZ_MASS
                    Dstr1_masd_CV[0] = (DP4_CV + PIS1P4_CV).M() - DP4_CV.M() + PDG_DZ_MASS

                    Dstr2_masd[0] = (EP4 + PIS2P4).M() - EP4.M() + PDG_DZ_MASS
                    Dstr2_masd_CV[0] = (EP4_CV + PIS2P4_CV).M() - EP4_CV.M() + PDG_DZ_MASS
                    # assigning info about the D E PIS
                    BV = TVector3(chain.B_DecayVtxX[cand], chain.B_DecayVtxY[cand], chain.B_DecayVtxZ[cand])
                    BVE = TVector3(sqrt(chain.B_DecayVtxXE[cand]), sqrt(chain.B_DecayVtxYE[cand]), sqrt(chain.B_DecayVtxZE[cand]))
                    BP4.SetXYZM( chain.B_px[cand], chain.B_py[cand], chain.B_pz[cand], chain.B_mass[cand])
                    BP3 = BP4.Vect()
                    B_pt[0]= BP4.Pt()
                    B_eta[0]= BP4.Eta()
                    B_vtxprob[0] = chain.B_Prob[cand]
                    B_mass[0]= chain.B_mass[cand]
                    B_masd[0]= chain.B_masd[cand]
                    B_masdd[0] = chain.B_masd[cand] - Dstr2_masd_CV[0] - Dstr1_masd_CV[0] + 2 * PDG_DSTR
                    B_pvdistsignif2[0] = DetachSignificance2( BV - PV, BVE, PVE)
                    B_pvdistsignif3[0] = DetachSignificance3( BV - PV, BVE, PVE)
                    B_pvdist[0] = (BV - PV).Mag()
                    B_pvcos2[0] = DirectionCos2 ( BV - PV, BP3)
                    B_pvcos3[0] = DirectionCos3 ( BV - PV, BP3)
                    if B_vtxprob[0] < 0.05 : continue

                    # assigning info about D to D E PIS
                    D_Bdistsignif2[0] = DetachSignificance2( DV - BV, DVE, BVE)
                    D_Bdistsignif3[0] = DetachSignificance3( DV - BV, DVE, BVE)
                    D_Bcos2[0] = DirectionCos2 ( DV - BV, DP3)
                    D_Bcos3[0] = DirectionCos3 ( DV - BV, DP3)

                    # assigning info about E to D E PIS
                    E_Bdistsignif2[0] = DetachSignificance2( EV - BV, EVE, BVE)
                    E_Bdistsignif3[0] = DetachSignificance3( EV - BV, EVE, BVE)
                    E_Bcos2[0] = DirectionCos2 ( EV - BV, EP3)
                    E_Bcos3[0] = DirectionCos3 ( EV - BV, EP3)

                    #cut on combinatorial effects in mischange pions and kaons - unnecessary
                    #if (EP4+PIS1P4).M() < 2.1 : continue

                    SAMEEVENT[0] = 0;
                    if (BBB > -1):
                        SAMEEVENT[0] = 1
                    #
                    BBB = 1; ## the next candidate is not the 1st in event
                    #
                    # fill _var branches in the tree if passes preselection cuts
                    mytree.Fill()
                #
                BBB = -1 ## when loop over candidates in event is finished, set this to -1, so the next candidate has SAMEEVENT=0
                if (evt % 20001 == 0):
                    print('Running... Now in event %i / %i ' %(evt, nevt))


            print('Total entries stored in MyTree: %s' %(mytree.GetEntries()))

            fileOUT.Write()
