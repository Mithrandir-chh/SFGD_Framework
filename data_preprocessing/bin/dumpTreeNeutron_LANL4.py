# created by Guang Yang Jan 2019
#!/usr/bin/env python

import sys
import os.path
import os
import ROOT
import random
import math
from optparse import OptionParser
from array import array

minEndInCM = [-12, -12, -50 ]
maxEndInCM = [12, 12, 50]

def getCube( start, end, posx, posy, posz ):

    pz= int(math.modf(posz)[1])
    py= int(math.modf(posy)[1])
    px= int(math.modf(posx)[1])

    res = []
    res.append(pz - start[2])
    res.append(py - start[1])
    res.append(px - start[0])
    res.append(res[0]*1000000 + res[1]*1000 + res[2])
    
    return res


def loop( events, tgeo, tout, nfiles, okruns ):

    flightDis = 90

    if len(okruns) == 0:
        print "There are no runs in this TTree...skipping!"
        return

    print "Inside event loop with %d files and first run %d" % (nfiles, okruns[0])

    event = ROOT.TG4Event()
    events.SetBranchAddress("Event",ROOT.AddressOf(event))
    print "Set branch address"

    freq = 3;
    current_time = -1;
    N = events.GetEntries()/10
    evt_per_file = N/nfiles
    if N % nfiles:
        print "Files don't all have the same number of events!!!"
        print "\n\n\n\n\n\n\n\n\n\n\n"

    print "Starting loop over %d entries" % N
    for ient in range(N):
        print '--------------------------------------------------------'
        print '--------------------------------------------------------'
        print '------------- A new event starts -----------------------'
        print '--------------------------------------------------------'
        print '--------------------------------------------------------'
        if ient % 100 == 0:
            print "Event %d of %d..." % (ient,N)
        events.GetEntry(ient)

        aRann = random.uniform(0,freq);
        if ient % int(347/freq) == 0:
            current_time += 1

        for ivtx,vertex in enumerate(event.Primaries):

            ## initialize output variables
            fileidx = ient/evt_per_file
            t_ifileNo[0] = okruns[fileidx]
            t_ievt[0] = ient%evt_per_file;
            t_vtx[0]=0.0; t_vtx[1]=0.0; t_vtx[2]=0.0;
	    for inii in range(1000):
	        t_neutronHitX[inii]=0.;
                t_neutronHitY[inii]=0.;
                t_neutronHitZ[inii]=0.;
                t_neutronHitT[inii]=0.;
                t_neutronHitE[inii]=0.;
                t_neutronHitS[inii]=0.;
                t_initialE[inii]=0.;
                t_daughterPx[inii]=0.;
                t_daughterPy[inii]=0.;
                t_daughterPz[inii]=0.;
                t_daughterPDG[inii]=0.;
                t_neutronPx[inii]=0.;
                t_neutronPy[inii]=0.;
                t_neutronPz[inii]=0.;
                t_node[inii]=-1;
                t_delay[inii]=-1;
                t_t0[inii]=-1;
                t_gamma_time[inii]=-1;
                t_neutron_time[inii]=-1;
                t_neutronId[inii]=0.;
                t_daughterParentId[inii]=0.;
                t_neutronParentId[inii]=0.;
                t_neutronParentPDG[inii]=0.;
		t_neutronTrackId[inii]=0.;
                t_neutronTrackPDG[inii]=0.;
                t_neutronTrackAng[inii]=0.;
		nTrackId[inii]=0;
		nParentId[inii]=0;
		nTrackEnergy[inii]=0;
		nTrackAng[inii]=0;
		nTrackPdg[inii]=0;
		nParentPdg[inii]=0;
		t_id[inii] = -1;
		t_length[inii] = 0.;
		if inii<100:
		    nelastic[inii] = 0;
            ## done

            # now ID numucc
            reaction=vertex.Reaction        

            # set the vertex location for output
            for i in range(3): 
                t_vtx[i] = vertex.Position[i] / 10. # cm
		#print t_vtx[i]

            nfsp = 0
            for ipart,particle in enumerate(vertex.Particles):
                e = particle.Momentum[3]
                p = (particle.Momentum[0]**2 + particle.Momentum[1]**2 + particle.Momentum[2]**2)**0.5
                #print "fs pdg list ",particle.PDGCode
                t_fsPdg[nfsp] = particle.PDGCode
                t_fsPx[nfsp] = particle.Momentum[0]
                t_fsPy[nfsp] = particle.Momentum[1]
                t_fsPz[nfsp] = particle.Momentum[2]
                t_fsE[nfsp] = e
                nfsp += 1

            t_nFS[0] = nfsp
	    # Get neutron trajectory information PDG neutron 2112 proton 2212
	    nTraj = []
            tag1=0;

	    incre1=0;
            for ctraj in event.Trajectories:
                if event.Trajectories[ctraj.TrackId].PDGCode != -1:

		    if ctraj.TrackId == 0:
		    	for cpoint in ctraj.Points:
			    if incre1 < 99:
			    	pointProcess[incre1]=cpoint.Process
                                pointSubProcess[incre1] = cpoint.Subprocess
			    	pointPositionX[incre1]=cpoint.Position.X()
			    	pointPositionY[incre1]=cpoint.Position.Y()
			    	pointPositionZ[incre1]=cpoint.Position.Z()
			    	incre1=incre1+1;

                    t_initialE[tag1] = event.Trajectories[ctraj.TrackId].InitialMomentum.Energy()
                    t_neutronPx[tag1] = event.Trajectories[ctraj.TrackId].InitialMomentum.Px()
                    t_neutronPy[tag1] = event.Trajectories[ctraj.TrackId].InitialMomentum.Py()
                    t_neutronPz[tag1] = event.Trajectories[ctraj.TrackId].InitialMomentum.Pz()
                    t_neutronId[tag1] = ctraj.TrackId
                    #print "a FS neutron with kinetic energy ", event.Trajectories[ctraj.TrackId].InitialMomentum.Energy()
                    tag1=tag1+1
                    #print ctraj.ParentId
                    #print ctraj.TrackId
                    #print ctraj.Name

            neutronMap = dict()
            tag2=0;
	    iid = 0;
            for ctraj in event.Trajectories:

		if event.Trajectories[ctraj.ParentId].PDGCode != -1:
                    t_daughterPx[tag2] = event.Trajectories[ctraj.TrackId].InitialMomentum.Px()
                    t_daughterPy[tag2] = event.Trajectories[ctraj.TrackId].InitialMomentum.Py()
                    t_daughterPz[tag2] = event.Trajectories[ctraj.TrackId].InitialMomentum.Pz()
                    t_daughterPDG[tag2] = event.Trajectories[ctraj.TrackId].PDGCode
                    t_daughterParentId[tag2] = ctraj.ParentId
                    #parentDaughterMap[ctraj.TrackId] = ctraj.ParentId
                    tag2=tag2+1
		    #print "a neutron daughter hit"
		    #print ctraj.ParentId
		    #print ctraj.TrackId
		    #print ctraj.Name
                    #print ctraj.InitialMomentum.Px(),' ',event.Trajectories[ctraj.TrackId].InitialMomentum.Px()
		    nTraj.append(ctraj.TrackId)
                    for ctraj4 in event.Trajectories:
                        if ctraj4.TrackId == ctraj.ParentId:
                            neutronMap[ctraj.TrackId] = ctraj4.TrackId                    
		    if iid< 100 :
		    	nTrackId[iid]     = ctraj.TrackId  
			nParentId[iid]    = ctraj.ParentId
			nTrackEnergy[iid] = ctraj.InitialMomentum.Energy()
			nTrackAng[iid]    = ctraj.InitialMomentum.Pz()/math.sqrt(pow(ctraj.InitialMomentum.Px(),2)+pow(ctraj.InitialMomentum.Px(),2)+pow(ctraj.InitialMomentum.Pz(),2))
			nTrackPdg[iid]	  = event.Trajectories[ctraj.TrackId].PDGCode
			nParentPdg[iid]	  = event.Trajectories[ctraj.ParentId].PDGCode
			iid += 1;

                #if ctraj.PDGCode == 2112:
		    #print ctraj.Name
                    #nTraj.append(ctraj.TrackId)

	    #print "----------------- ",ient

	    # neutron hit number looper that will be used soon
	    iN = 0;            

	    # hits in 3DST 

            hits = []
            for key3 in event.SegmentDetectors:
                #if "volCube" in key3.first or "volcube" in key3.first or "" in key3.first:
		if "volCube" in key3.first:
                    hits += key3.second

	    recordTrackId = 100*[-1]

            for hit in hits:
                hStart = ROOT.TVector3( hit.Start[0]/10., hit.Start[1]/10., hit.Start[2]/10. )
                hStop = ROOT.TVector3( hit.Stop[0]/10., hit.Stop[1]/10., hit.Stop[2]/10. )

		#print hit.Contrib[0]
		#Neutron PDG code 2112, neutron should not give energy deposit itself
	 	for iloop in nTraj:
		    #print "checking"
		    #print hit.Contrib[0],iloop
		    if hit.Contrib[0] == iloop and iN < 999:

                        for innLoop in range(100):
                            if hit.EnergyDeposit > innLoop*0.05 and recordTrackId[innLoop] == -1:
                                recordTrackId[innLoop] = event.Trajectories[iloop].TrackId

			t_id[iN] = iloop
			t_length[iN] = (hStop-hStart).Mag()
			#print hStart.x()
			#print event.Trajectories[iloop].Name
			t_neutronHitX[iN] = hStart.x()
                        t_neutronHitY[iN] = hStart.y()
                        t_neutronHitZ[iN] = hStart.z()
                        #nodeInfo   = str(tgeo.FindNode( hStart.x(),hStart.y(),hStart.z() ).GetName())
                        #t_node[iN] = int(filter(str.isdigit, nodeInfo))
                        tgeo.FindNode( hStart.x(),hStart.y(),hStart.z() )
                        cubeInfo          = getCube(minEndInCM, maxEndInCM, (hStart.x()+hStop.x())/2., (hStart.y()+hStop.y())/2., (hStart.z()+hStop.z())/2.)
                        #print 'location of the hit ',cubeInfo[0],' ',cubeInfo[1],' ',cubeInfo[2],' ',cubeInfo[3]                       

                        #t_node[iN]        = cubeInfo[3]

                        t_node[iN] = tgeo.GetCurrentNodeId()
                    
                        #print tgeo.FindNode( hStart.x(),hStart.y(),hStart.z()).GetName(), tgeo.GetCurrentNodeId()                                            

                        t_neutronHitT[iN] = hit.Start.T()
                        t_neutronHitE[iN] = hit.EnergyDeposit
                        #print 'energy deposit ',t_neutronHitE[iN]
                        #print 'node is ',t_node.GetName()
                        t_neutronHitS[iN] = event.Trajectories[iloop].PDGCode
                        #if(event.Trajectories[iloop].PDGCode>10000):
                            #print 'found a bad PDG *********************** ', event.Trajectories[iloop].Name

                        # time structure: gtrig 10 us on-going
                        # spill start for 650 us
                        # neutron hit time in 1.8 us
                        # assume 5 x 1.8us there is one neutron

                        #print current_time, ' ',ient,' ', (ient - current_time * int(347/freq))
                        t_t0[iN] = int(1800 * int((ient - current_time * int(347/freq)) * freq + aRann));

                        t_spill[iN] = current_time; 

                        t_gamma_time[iN] = t_t0[iN] + flightDis/0.3  ;
                        # double d = 90;
                        # double e = i/10.;
                        # double m = 939;
                        # double c = 0.3;

                        # double gt = d/c;
                        # double gamm  = e/m + 1.;
                        # double v = sqrt(-( (1/gamm)*(1/gamm) * c * c - c * c));
                        # double t = d / v;

                        # TOF = 72.3 * L / sqrt(E) ; E in MeV, TOF in ns, L in m
                        
                        #tof = 72.3 * flightDis / math.sqrt(t_initialE[0]) ;

                        d = 90;
                        e = t_initialE[0];
                        m = 939;
                        c = 0.3;

                        gt = d/c;
                        gamm  = e/m + 1.;
                        v = math.sqrt(-( (1/gamm)*(1/gamm) * c * c - c * c));
                        tof = d / v;

                        t_neutron_time[iN] = t_t0[iN] + tof;

                        #print 't_t0 ', t_t0[iN], ' t_spill[iN] ',t_spill[iN], ' t_gamma_time ',t_gamma_time[iN], ' t_neutron_time ', t_neutron_time[iN]

                        #print 'initial neutron energy ',t_initialE[0], 'neutron tof is ',tof, ' delay number ',int(tof/1800)
                        t_delay[iN] = int(tof/1800);
                        if (neutronMap.get(iloop) > -1000000):
                            #t_neutronParentId[iN] = neutronMap.get(iloop)
                            #t_neutronParentPDG[iN] = event.Trajectories[neutronMap.get(iloop)].PDGCode
                            t_neutronParentId[iN] = event.Trajectories[iloop].ParentId;
                            t_neutronParentPDG[iN] = event.Trajectories[event.Trajectories[iloop].ParentId].PDGCode;
			    t_neutronTrackId[iN] = event.Trajectories[iloop].TrackId;
			    t_neutronTrackPDG[iN] = event.Trajectories[iloop].PDGCode; 
			    t_neutronTrackAng[iN] = event.Trajectories[iloop].InitialMomentum.Pz()/math.sqrt(pow(event.Trajectories[iloop].InitialMomentum.Px(),2)+pow(event.Trajectories[iloop].InitialMomentum.Px(),2)+pow(event.Trajectories[iloop].InitialMomentum.Pz(),2));
 
                        iN += 1;

            ## end the hit loop here

            #for inloop in range(100):
		#if recordTrackId[inloop] > -1:
		#    Nelastic[inloop] = 0
                #for ctraj in event.Trajectories:

                    #if event.Trajectories[ctraj.ParentId].PDGCode != -1 and recordTrackId[inloop] == ctraj.TrackId:
                    #    t_daughterPx[tag2] = event.Trajectories[ctraj.TrackId].InitialMomentum.Px()
                    #    t_daughterPy[tag2] = event.Trajectories[ctraj.TrackId].InitialMomentum.Py()
                    #    t_daughterPz[tag2] = event.Trajectories[ctraj.TrackId].InitialMomentum.Pz()
                    #    t_daughterPDG[tag2] = event.Trajectories[ctraj.TrackId].PDGCode
                    #    t_daughterParentId[tag2] = ctraj.ParentId
                        #parentDaughterMap[ctraj.TrackId] = ctraj.ParentId
                    #    tag2=tag2+1
                        #print " ::::::::: the first hit location, trackId, parentId, trackName, parentName, energy, angle, pdg, parentPdg : ",inloop," ",ctraj.TrackId," ",ctraj.ParentId," ",ctraj.Name," ",event.Trajectories[ctraj.ParentId].Name," ",ctraj.InitialMomentum.Energy()," ",ctraj.InitialMomentum.Pz()/math.sqrt(pow(ctraj.InitialMomentum.Px(),2)+pow(ctraj.InitialMomentum.Px(),2)+pow(ctraj.InitialMomentum.Pz(),2))," ",event.Trajectories[ctraj.TrackId].PDGCode," ",event.Trajectories[ctraj.ParentId].PDGCode

                        #saveID = 20* [-1]
                        #if recordTrackId[inloop] == 0:
                        #    Nelastic[inloop] = -1
			#elif recordTrackId[inloop] == 1:
                        #    Nelastic[inloop] = 0
                        #elif ctraj.ParentId == 0 :
                        #    Nelastic[inloop] = recordTrackId[inloop]-1
                        #elif ctraj.ParentId == 1 :
                        #    Nelastic[inloop] = recordTrackId[inloop]-1
                        #else:
                        #    tem4p = ctraj.TrackId
                        #    for kkk in range(20):
                        #        if event.Trajectories[tem4p].ParentId == -1 or event.Trajectories[tem4p].TrackId == 0:
                        #            break;
                        #        saveID[kkk] = event.Trajectories[tem4p].ParentId
                        #        tem4p = event.Trajectories[tem4p].ParentId;

                        #for ctrajj in event.Trajectories:
                        #    for mmm in range(20):
                        #        if saveID[mmm] == -1:
                        #            break;
                        #        else:
                        #            if ctrajj.ParentId == event.Trajectories[saveID[mmm]].TrackId and ctrajj.TrackId < recordTrackId[inloop]:
                        #                Nelastic[inloop] += 1
			#if Nelastic[inloop]> 0 :
                            #print " >>>>>>>>>>>>>>>> testing result of number of scattering : ", Nelastic[inloop] , " more checks of saveID[0-4] : ",saveID[0]," ",saveID[1]," ",saveID[2]," ",saveID[3] 
	    #for inloop in range(100):
	    	#if Nelastic[inloop] > 0:
		    #print " >>>>>>>>>>>>>>>> testing result of number of scattering : ", Nelastic[inloop], "  at  ",inloop

            print ' ------------------ Time for track details ------------------------'
            print ' ------------------------------------------------------------------'
            for ctraj in event.Trajectories:
                if event.Trajectories[ctraj.ParentId].PDGCode != -1:
                    print " ::: trackId, parentId, trackName, parentName, energy, angle, pdg, parentPdg : ",ctraj.TrackId," ",ctraj.ParentId," ",ctraj.Name," ",event.Trajectories[ctraj.ParentId].Name," ",ctraj.InitialMomentum.Energy()," ",ctraj.InitialMomentum.Pz()/math.sqrt(pow(ctraj.InitialMomentum.Px(),2)+pow(ctraj.InitialMomentum.Px(),2)+pow(ctraj.InitialMomentum.Pz(),2))," ",event.Trajectories[ctraj.TrackId].PDGCode," ",event.Trajectories[ctraj.ParentId].PDGCode

                    print "track points details.. "
                    for cpoint in ctraj.Points:
                        print "point process, subprocess and position : ",cpoint.Process, " ",cpoint.Subprocess," ",  cpoint.Position.X()," ",cpoint.Position.Y()," ",cpoint.Position.Z()
                        print "point position.energy, point.momentum.mag2 and point.momentum.xyz : ", cpoint.Position.Energy()," ", cpoint.Momentum.Mag2()," ",cpoint.Momentum.Px()," ",cpoint.Momentum.Py()," ",cpoint.Momentum.Pz()


            print ' ------------------------------------------------------------------'
            print ' ------------------------------------------------------------------'
            print ' ----------------------- # elastic calculation --------------------'
            print ' ------------------------------------------------------------------'
            print ' ------------------------------------------------------------------'


            #nelastic = 100* [0.]
            #for hit in hits:

            for inloop in range(50):
                if nelastic[inloop] != 0:
                    continue;
                zbarr = 240 + inloop*10 
                recordedLocation = 100*[-2.]
                for ctraj in event.Trajectories:

                    existing = False;
                    for trajPoint in event.Trajectories[ctraj.TrackId].Points:
                        if trajPoint.Position.Z() >= zbarr:
                            break;
                        elif trajPoint.Process == 4:
                            for inner in range(int(nelastic[inloop])):
                                if recordedLocation[inner] == trajPoint.Position.Z():
                                    existing = True;
                            if existing == True:
                                continue;
                            recordedLocation[int(nelastic[inloop])] = trajPoint.Position.Z();
                            nelastic[inloop] += 1;

                    #print "with energy at ",inloop," , elastic scattering number : ",nelastic[inloop], "  from track : ",recordTrackId[inloop]

            tout.Fill()

if __name__ == "__main__":

    ROOT.gROOT.SetBatch(1)

    parser = OptionParser()
    parser.add_option('--outfile', help='Output file name', default="out.root")
    parser.add_option('--topdir', help='Input file top directory', default="")
    parser.add_option('--first_run', type=int, help='First run number', default=0)
    parser.add_option('--last_run', type=int, help='Last run number', default=0)
    parser.add_option('--rhc', action='store_true', help='Reverse horn current', default=False)
    parser.add_option('--grid', action='store_true', help='grid mode', default=False)

    (args, dummy) = parser.parse_args()

    # make an output ntuple
    fout = ROOT.TFile( args.outfile, "RECREATE" )
    tout = ROOT.TTree( "tree","tree" )
    t_ifileNo = array('i',[0])
    tout.Branch('ifileNo',t_ifileNo,'ifileNo/I')
    t_ievt = array('i',[0])
    tout.Branch('ievt',t_ievt,'ievt/I')

    t_vtx = array('f',3*[0.0])
    tout.Branch('vtx',t_vtx,'vtx[3]/F')

    t_nFS = array('i',[0])
    tout.Branch('nFS',t_nFS,'nFS/I')
    t_fsPdg = array('i',100*[0])
    tout.Branch('fsPdg',t_fsPdg,'fsPdg[nFS]/I')
    t_fsPx = array('f',100*[0.])
    tout.Branch('fsPx',t_fsPx,'fsPx[nFS]/F')
    t_fsPy = array('f',100*[0.])
    tout.Branch('fsPy',t_fsPy,'fsPy[nFS]/F')
    t_fsPz = array('f',100*[0.])
    tout.Branch('fsPz',t_fsPz,'fsPz[nFS]/F')
    t_fsE = array('f',100*[0.])
    tout.Branch('fsE',t_fsE,'fsE[nFS]/F')

    t_neutronHitX = array('f' ,1000*[0.])
    tout.Branch('neutronHitX',t_neutronHitX,'neutronHitX[1000]/F')
    t_neutronHitY = array('f' ,1000*[0.])
    tout.Branch('neutronHitY',t_neutronHitY,'neutronHitY[1000]/F')
    t_neutronHitZ = array('f' ,1000*[0.])
    tout.Branch('neutronHitZ',t_neutronHitZ,'neutronHitZ[1000]/F')
    t_neutronHitT = array('f' ,1000*[0.])
    tout.Branch('neutronHitT',t_neutronHitT,'neutronHitT[1000]/F')
    t_neutronHitE = array('f' ,1000*[0.])
    tout.Branch('neutronHitE',t_neutronHitE,'neutronHitE[1000]/F')
    t_neutronHitS = array('f' ,1000*[0.])
    tout.Branch('neutronHitPDG',t_neutronHitS,'neutronHitPDG[1000]/F')
    t_initialE = array('f' ,1000*[0.])
    tout.Branch('initialE',t_initialE,'initialE[1000]/F')
    t_daughterPx = array('f' ,1000*[0.])
    tout.Branch('daughterPx',t_daughterPx,'daughterPx[1000]/F')
    t_daughterPy = array('f' ,1000*[0.])
    tout.Branch('daughterPy',t_daughterPy,'daughterPy[1000]/F')
    t_daughterPz = array('f' ,1000*[0.])
    tout.Branch('daughterPz',t_daughterPz,'daughterPz[1000]/F')
    t_daughterPDG = array('f' ,1000*[0.])
    tout.Branch('daughterPDG',t_daughterPDG,'daughterPDG[1000]/F')
    t_node = array('i' ,1000*[0])
    tout.Branch('node',t_node,'node[1000]/I')
    t_neutronPx = array('f' ,1000*[0.])
    tout.Branch('neutronPx',t_neutronPx,'neutronPx[1000]/F')
    t_neutronPy = array('f' ,1000*[0.])
    tout.Branch('neutronPy',t_neutronPy,'neutronPy[1000]/F')
    t_neutronPz = array('f' ,1000*[0.])
    tout.Branch('neutronPz',t_neutronPz,'neutronPz[1000]/F')
    t_neutronId = array('f' ,1000*[0.])
    tout.Branch('neutronId',t_neutronId,'neutronId[1000]/F')
    t_daughterParentId = array('f' ,1000*[0.])
    tout.Branch('daughterParentId',t_daughterParentId,'daughterParentId[1000]/F')
    t_delay = array('i' ,1000*[0])
    tout.Branch('delay',t_delay,'delay[1000]/I')
    t_t0 = array('f', 1000*[0])
    tout.Branch('t0',t_t0, 't_t0[1000]/F')
    t_spill = array('f', 1000*[0])
    tout.Branch('spill',t_spill, 't_spill[1000]/F')
    t_gamma_time = array('f', 1000*[0])
    tout.Branch('gamma_time',t_gamma_time,'t_gamma_time[1000]/F')
    t_neutron_time = array('f', 1000*[0])
    tout.Branch('neutron_time', t_neutron_time, 't_neutron_time')

    t_neutronParentId = array('f' ,1000*[0.])
    tout.Branch('neutronParentId',t_neutronParentId,'neutronParentId[1000]/F')
    t_neutronParentPDG = array('f' ,1000*[0.])
    tout.Branch('neutronParentPDG',t_neutronParentPDG,'neutronParentPDG[1000]/F')
    t_neutronTrackId = array('f' ,1000*[0.])
    tout.Branch('neutronTrackId',t_neutronTrackId,'neutronTrackId[1000]/F')
    t_neutronTrackPDG = array('f' ,1000*[0.])
    tout.Branch('neutronTrackPDG',t_neutronTrackPDG,'neutronTrackPDG[1000]/F')
    t_neutronTrackAng = array('f' ,1000*[0.])
    tout.Branch('neutronTrackAng',t_neutronTrackAng,'neutronTrackAng[1000]/F')

    t_id = array('i', 1000*[0])
    tout.Branch('id',t_id, 't_id[1000]/I')
    t_length = array('f', 1000*[0])
    tout.Branch('length',t_length, 't_length[1000]/F')

    nTrackId = array('i', 1000*[0])
    tout.Branch('nTrackId',nTrackId, 'nTrackId[1000]/I')
    nParentId = array('i', 1000*[0])
    tout.Branch('nParentId',nParentId, 'nParentId[1000]/I')
    nTrackEnergy = array('f', 1000*[0])
    tout.Branch('nTrackEnergy',nTrackEnergy, 'nTrackEnergy[1000]/F')
    nTrackAng = array('f', 1000*[0])
    tout.Branch('nTrackAng',nTrackAng, 'nTrackAng[1000]/F')
    nTrackPdg = array('i', 1000*[0])
    tout.Branch('nTrackPdg',nTrackPdg, 'nTrackPdg[1000]/I')
    nParentPdg = array('i', 1000*[0])
    tout.Branch('nParentPdg',nParentPdg, 'nParentPdg[1000]/I')

    pointPositionX = array('f', 100*[0])
    tout.Branch('pointPositionX',pointPositionX, 'pointPositionX[100]/F')
    pointPositionY = array('f', 100*[0])
    tout.Branch('pointPositionY',pointPositionY, 'pointPositionY[100]/F')
    pointPositionZ = array('f', 100*[0])
    tout.Branch('pointPositionZ',pointPositionZ, 'pointPositionZ[100]/F')
    pointProcess = array('i', 100*[0])
    tout.Branch('pointProcess',pointProcess, 'pointProcess[100]/I')
    pointSubProcess = array('i', 100*[0])
    tout.Branch('pointSubProcess',pointSubProcess, 'pointSubProcess[100]/I')

    nelastic = array('f', 100*[0.])
    tout.Branch('nElastic',nelastic,'nElastic[100]/F')

    loaded = False
    #if os.path.exists("EDepSimEvents/EDepSimEvents.so"):
    #    print "Found EDepSimEvents.so"
    #    ROOT.gSystem.Load("EDepSimEvents/EDepSimEvents.so")
    #    loaded = True

    tgeo = None

    events = ROOT.TChain( "EDepSimEvents", "main event tree" )

    neutrino = "neutrino"
    if args.rhc:
        neutrino = "antineutrino"

    print "Building TChains for runs %d-%d..." % (args.first_run, args.last_run)
    nfiles = 0
    okruns = []
    for run in range( args.first_run, args.last_run+1 ):
        if args.grid == False:
            fname = "%s/edep.%d.root" % (args.topdir,run)
        else:
            fname = "%s/%02d/full3DST.%s.%d.edepsim.root" % (args.topdir, run/1000, neutrino, run)
        print fname

        # see if it is an OK file
        if not os.access( fname, os.R_OK ):
            print "Can't access file: %s" % fname
            continue
        tf = ROOT.TFile( fname )
        if tf.TestBit(ROOT.TFile.kRecovered): # problem with file
            print "File is crap: %s" % fname
            continue
        nfiles += 1
        okruns.append( run )

        if not loaded:
            loaded = True
            tf.MakeProject("EDepSimEvents","*","RECREATE++")

        # add it to the tchain
        events.Add( fname )

        if tgeo is None: # first OK file, get geometry
            tgeo = tf.Get("EDepSimGeometry")
        tf.Close() # done with this one

    print "OK runs: ", sorted(okruns)
    print "got %d events in %d files = %1.1f events per file" % (events.GetEntries(), nfiles, 1.0*events.GetEntries()/nfiles)
    loop( events, tgeo, tout, nfiles, sorted(okruns) )

    fout.cd()
    tout.Write()




