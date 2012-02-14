from numpy import *

def findtheinitial(oplfile, pickupfile):
    if oplfile[:4].isdigit():
	outputfile = "InitialPosition" + oplfile[:4] +".txt"
    else:
	outputfile = "InitialPosition.txt"
    opl=genfromtxt(oplfile)
    pickup=genfromtxt(pickupfile).astype('int64')
    init = repeat(0,4)
    for i in r_[:pickup.size]:
        meet, = (opl[:,0]==pickup[i]).nonzero()
        if meet.size > 0:
            init = vstack([init, opl[meet[0], :4]])
    init = init[1:,:]
    savetxt(outputfile, init, fmt='%d\t%d\t%10.5f\t%10.5f')        

def findtheterminal(oplfile, pickupfile):
    if oplfile[:4].isdigit():
	outputfile = "TerminalPosition" + oplfile[:4] +".txt"
    else:
	outputfile = "TerminalPosition.txt"
    opl=genfromtxt(oplfile)
    pickup=genfromtxt(pickupfile).astype('int64')
    term = repeat(0,4)
    for i in r_[:pickup.size]:
        meet, = (opl[:,0]==pickup[i]).nonzero()
        if meet.size > 0:
            term = vstack([term, opl[meet[-1],:4]])
    term = term[1:,:]
    savetxt(outputfile, term, fmt='%d\t%d\t%10.5f\t%10.5f')

def findthebothends(oplfile, pickupfile):
    if oplfile[:4].isdigit():
	outputfile = "BothEnds" + oplfile[:4] +".txt"
    else:
	outputfile = "TerminalPosition.txt"
    opl=genfromtxt(oplfile)
    pickup=genfromtxt(pickupfile).astype('int64')
    ends = repeat(0,7)
    for i in r_[:pickup.size]:
        meet, = (opl[:,0]==pickup[i]).nonzero()
        if meet.size > 0:
            ends = vstack([ends, r_[opl[meet[0],:4], opl[meet[-1],1:4]]])
            #        ends[i,:4] = opl[meet[0],:4]
            #        ends[i,4:8] = opl[meet[-1],1:4]
    ends = ends[1:,:]
    savetxt(outputfile, ends, fmt='%d\t%d\t%10.5f\t%10.5f\t%d\t%10.5f\t%10.5f')

def removeoutsidetracks(oplfile):
    "Use results table just tagged OnCell in ImageJ as an input file"
    if oplfile.find("OnCell"):
	outputfile = oplfile.replace("OnCell", "cell")
    else:
	outputfile = oplfile[:-11]+"_cell.txt"
    opl=genfromtxt(oplfile,skiprows=1)
    idoutside=unique(opl[opl[:,5]==0,1])
    abolish=asarray([opl[i,1] in idoutside for i in r_[:opl.shape[0]]])
    oplinside=opl[abolish==False,1:5]
    savetxt(outputfile, oplinside, fmt='%d\t%d\t%10.5f\t%10.5f')

def removeoutsidetrackswithdensity(oplfile):
    "Use results table just tagged OnCell in ImageJ as an input file"
    if oplfile.find("OnCell"):
	outputfile = oplfile.replace("OnCell", "cell")
    else:
	outputfile = oplfile[:-11]+"_cell.txt"
    opl=genfromtxt(inoplfilewtag,skiprows=1)
    idoutside=unique(opl[opl[:,6]==0,1])
    abolish=asarray([opl[i,1] in idoutside for i in r_[:opl.shape[0]]])
    oplinside=opl[abolish==False,1:6]
    savetxt(outputfile, oplinside, fmt='%d\t%d\t%10.5f\t%10.5f\t%d')

def pickupopltail(oplfile, pickupfile):
    outputfile = oplfile[:-4]+"_TerTrk.txt"
    minlen=22
    opl=genfromtxt(oplfile)
    pickup=genfromtxt(pickupfile).astype('int64')
    trk = repeat(0,4)
    for i in r_[:pickup.size]:
	    #	print "i=",i,pickup[i]
	    #        print (opl[:,0]==pickup[i])
        meet, = (opl[:,0]==pickup[i]).nonzero()
        if meet.size >= minlen:
		#	    print meet[0],meet[-1],opl[meet[0]],opl[meet[-1]]
            trk = vstack([trk, r_[opl[meet[-minlen]:meet[-1]+1,:4]]])
    trk = trk[1:,:]
    savetxt(outputfile, trk, fmt='%d\t%d\t%10.5f\t%10.5f')

def makeoplreference(oplfile):
    outputfile = oplfile[:-4]+"_ref.txt"	
    opl=genfromtxt(oplfile)
    ids = unique(opl[:,0])
    ends = repeat(0,7)
    for id in ids:
	meet, = (opl[:,0]==id).nonzero()
        if meet.size > 0:
	    ends = vstack([ends, r_[opl[meet[0],:4], opl[meet[-1],1:4]]])
    ends = ends[1:,:]
    savetxt(outputfile, ends, fmt='%d\t%d\t%10.5f\t%10.5f\t%d\t%10.5f\t%10.5f')

# def fixdiscontinuouspoint(oplfile,outputfile):
#     thr = 5
#     opl=genfromtxt(oplfile)
#     def sqd(v1,v2):
# 	return (v1[2]-v2[2])**2 + (v1[3]-v2[3])**2
    
#     sqthr=thr**2
#     for i in r_[1:opl.shape[0]-1]:
# 	if (opl[i-1,0]==opl[i+1,0]):
# 	    if (sqd(opl[i,:],opl[i-1,:])>sqthr) and (sqd(opl[i,:],opl[i+1,:])>sqthr) and (sqd(opl[i-1,:],opl[i+1,:])<4*sqthr):
# 		opl[i,2]=(opl[i-1,2]+opl[i+1,2])/2.0
# 		opl[i,3]=(opl[i-1,3]+opl[i+1,3])/2.0
#     savetxt(outputfile, opl, fmt='%d\t%d\t%10.5f\t%10.5f')
	
