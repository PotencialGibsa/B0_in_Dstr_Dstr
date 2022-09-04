from ROOT import *
import datetime
###################
#SimpleFiles_DD_v4.root_16_04_connected.root
#'Simplefiles_DD_v4_2018_connected.root'
#'Simplefiles_SiEl_Dst_D0_2017par_v4_connected.root'
time = datetime.datetime.now().timetuple()
time_st = str(time[3])+'_'+str(time[4])+'_'+str(time[5])+'_'

ch = TChain('mytree')
#ch.Add('Simplefiles_SiEl_Dst_D0_2017par_v4_connected.root')

ch.Add('/eos/home-d/dshmygol/DstDstv1/SimpleFile_DstDstv1_A.root')
ch.Add('/eos/home-d/dshmygol/DstDstv1/SimpleFile_DstDstv1_B.root')
ch.Add('/eos/home-d/dshmygol/DstDstv1/SimpleFile_DstDstv1_C.root')
ch.Add('/eos/home-d/dshmygol/DstDstv1/SimpleFile_DstDstv1_D1.root')
ch.Add('/eos/home-d/dshmygol/DstDstv1/SimpleFile_DstDstv1_D2.root')
ch.Add('/eos/home-d/dshmygol/DstDstv1/SimpleFile_DstDstv1_D3.root')
ch.Add('/eos/home-d/dshmygol/DstDstv1/SimpleFile_DstDstv1_D4.root')
ch.Add('/eos/home-d/dshmygol/DstDstv1/SimpleFile_DstDstv1_D5.root')

nEvt = ch.GetEntries()
Bmass_intr = [5.0, 6.0]
print(nEvt)

hfile = TFile( 'hist-DstDstv1_'+time_st+'.root', 'RECREATE', 'Demo ROOT file with histograms' )
hDxDs_cc = TH2D('D* x D*_cc', 'D* x D*_cc', 100, 2.005, 2.024 , 100 , 2.005, 2.024)
hDxDs_c_anti_c = TH2D('D* x D*_c_anti_c', 'D* x D*_c_anti_c', 100, 1.855, 1.87 , 100 , 2.005, 2.024)
hBm_cc = TH1F("B_masd_cc", "B_masd_cc", 50, 4.0, 4.1)
hBm_c_anti_c = TH1F("B_masd_c_anti_c", "B_masd_c_anti_c", 50, 3.87, 4.5)
hBm = TH1F("B_masd", "B_masd", 60, Bmass_intr[0], Bmass_intr[1])
n = 0

mass_cmp = []
comb_n = 0

for evt in range(nEvt):
    if ch.GetEntry(evt) <= 0: break

    E_mass = ch.E_mass
    D_mass = ch.D_mass
    Dstr1_masd = ch.Dstr1_masd
    Dstr2_masd = ch.Dstr2_masd
    E_pt = ch.E_pt
    E_pvdistsignif2 = ch.E_pvdistsignif2
    E_pvcos2 = ch.E_pvcos2
    #B_masd = ch.B_masd - ch.Dstr1_masd - ch.E_mass + 1.86483 + 2.010
    B_masd = ch.B_masd
    PIS1_ch = ch.PIS1_chrg
    PIS2_ch = ch.PIS2_chrg
    KA2_pt = ch.KA2_pt
    PI2_pt = ch.PI2_pt
    KA2_ch = ch.KA2_chrg
    B_pvdistsignif2 = ch.B_pvdistsignif2

    if B_masd > Bmass_intr[1]: continue
    if B_masd < Bmass_intr[0]: continue
    if B_pvdistsignif2 < 3 : continue
    if ch.B_pvcos2 < 0.9 : continue
    if Dstr1_masd < 2.006: continue
    if Dstr1_masd > 2.020: continue
    if Dstr2_masd < 2.006: continue
    if Dstr2_masd > 2.020: continue

    n+=1
    print(n)
    #if n%10 == 0:
        #print(n, "nevt")

    #comb
    elem_mass_ar = [ch.KA1_pt , ch.KA2_pt , ch.PI1_pt, ch.PI2_pt, ch.PIS1_pt]
    comb_condition = False
    '''
    for iter in mass_cmp:
            intersection = list(set(elem_mass_ar) & set(iter))
            for intersection_iter in intersection:
                if elem_mass_ar.index(intersection_iter) != iter.index(intersection_iter):
                    comb_condition = True
                    continue
            if comb_condition:
                continue

        if comb_condition:
            comb_n += 1
            if comb_n%1000 == 0:
                print(comb_n, "comb_n")

            continue
        else:
            mass_cmp.append(sorted([ch.KA1_pt , ch.KA2_pt , ch.PI1_pt, ch.PI2_pt, ch.PIS_pt]))
    '''
        ##Fill
    '''
        if PIS1_ch * KA2_ch > 0:
            if E_mass <  1.869486813 - 0.0125963048779 * 3: continue #	if E_mass < 1.76: continue
            if E_mass > 1.869486813 + 0.0125963048779 * 3: continue
            hDxDs_cc.Fill(D_mass , Dstr1_masd) # or E_mass
            hBm_cc.Fill(B_masd)
            hBm.Fill(B_masd)

        if PIS_ch * KA2_ch < 0:
            if E_mass <  1.863182 - 0.011357*3: continue #	if E_mass < 1.76: continue
            if E_mass > 1.863182 + 0.011357*3: continue
            hDxDs_c_anti_c.Fill(D_mass , Dstr1_masd)
            hBm_c_anti_c.Fill(B_masd)
            hBm.Fill(B_masd)
    '''
        	#print(comb_n)
    hBm.Fill(B_masd)

canvas = TCanvas("canvas")
canvas.cd()

hDxDs_cc.SetTitle('D* & D0 cc')
hDxDs_cc.SetYTitle('M(D*), [GeV]')
hDxDs_cc.SetXTitle('M(D0), [GeV]')
hDxDs_cc.SetLabelOffset (0.005, axis = "Y")
hDxDs_cc.SetTitleOffset (1.35 , axis = "Y")
hDxDs_cc.SetTitleOffset (1.35 , axis = "X")
hDxDs_cc.Draw()

canvas.SaveAs("hDxDs_cc.png")

hBm_cc.SetTitle('B_masd cc')
hBm_cc.SetYTitle('Candidates')
hBm_cc.SetXTitle('M(B_masd), [GeV]')
hBm_cc.SetLabelOffset (0.005, axis = "Y")
hBm_cc.SetTitleOffset (1.35 , axis = "Y")
hBm_cc.SetTitleOffset (1.35 , axis = "X")
hBm_cc.Draw()

canvas.SaveAs("hBm_cc.png")

hDxDs_c_anti_c.SetTitle('D* & D0 c_anti_c')
hDxDs_c_anti_c.SetYTitle('M(D*), [GeV]')
hDxDs_c_anti_c.SetXTitle('M(D0), [GeV]')
hDxDs_c_anti_c.SetLabelOffset (0.005, axis = "Y")
hDxDs_c_anti_c.SetTitleOffset (1.35 , axis = "Y")
hDxDs_c_anti_c.SetTitleOffset (1.35 , axis = "X")
hDxDs_c_anti_c.Draw()

canvas.SaveAs("hDxDs_c_anti_c.png")

hBm_c_anti_c.SetTitle('B_masd c anti c')
hBm_c_anti_c.SetYTitle('Candidates')
hBm_c_anti_c.SetXTitle('M(B_masd), [GeV]')
hBm_c_anti_c.SetLabelOffset (0.005, axis = "Y")
hBm_c_anti_c.SetTitleOffset (1.35 , axis = "Y")
hBm_c_anti_c.SetTitleOffset (1.35 , axis = "X")
hBm_c_anti_c.Draw()

canvas.SaveAs("hBm_c_anti_c.png")

hBm.SetTitle('B_masd')
hBm.SetYTitle('Candidates')
hBm.SetXTitle('M(B_masd), [GeV]')
hBm.SetLabelOffset (0.005, axis = "Y")
hBm.SetTitleOffset (1.35 , axis = "Y")
hBm.SetTitleOffset (1.35 , axis = "X")
hBm.Draw()

canvas.SaveAs("hBm.png")

hfile.Write()
