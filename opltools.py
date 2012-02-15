from numpy import *
import os

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

def concatpolishdir(folderpath):
    """Concatenate Polished Spots files for reverse indexing"""
    xyzs = zeros([1,6])
    for file in os.listdir(folderpath):
        if file.startswith('PSpA'):
            frame=int(file[4:8])
            filepath=folderpath+'/'+file
            data=genfromtxt(filepath,skiprows=1)
            data=data[data[:,12]==1,:6]
            xyzs=row_stack([xyzs,data])
    return xyzs[1:,:]

# def reversemap(oplfile,polishdir):
#     outputfile = oplfile[:-4]+"_revmap.txt"
#     opl=genfromtxt(oplfile)
#     xyzs=concatpolishdir(polishdir)
#     tagged = zeros([1,4])
#     for j in r_[:xyzs.shape(0)]:
#         oid, X, Y, slice, cx, cy = xyzs[j,:]
#         target = ( (opl[:,1]==slice) and (opl[:,2]==cx) and (opl[:,3]==cy) )
#         sid = opl[target,0]
#         tagged = vstack([tagged, r_[cx,cy,slice,sid]])
#     tagged = tagged[1:,:]
#     savetxt("xyzs.txt", tagged, fmt='%10.5f\t%10.5f\t%d\t%d')
    
def backtrace(oplfile,polishdir):
    outputfile = oplfile[:-4]+"_backtraced.txt"
    opl=genfromtxt(oplfile)
    xyzs=concatpolishdir(polishdir)
    oplcolnum = opl.shape[1]
    if oplcolnum == 4:
        tagged = zeros([1,oplcolnum+1])
        for j in r_[:opl.shape[0]]:
            sid, slice, cx, cy = opl[j,:4]
            #            target = (xyzs[:,3]==slice)
            target = ( (xyzs[:,3]==slice) & (xyzs[:,4]==cx) & (xyzs[:,5]==cy) )
            if target.any():
                subid = xyzs[target,0]
            else:
                subid = -1
            tagged = vstack([tagged, r_[opl[j,:],subid]])
        tagged = tagged[1:,:]
        savetxt(outputfile, tagged, fmt='%d\t%d\t%10.5f\t%10.5f\t%d')

def backtraceref(oplfile,polishdir):
    outputfile = oplfile[:-4]+"_backtracedref.txt"
    opl=genfromtxt(oplfile)
    xyzs=concatpolishdir(polishdir)
    oplcolnum = opl.shape[1]
    if oplfile.endswith("_ref.txt") and oplcolnum == 7:
        tagged = zeros([1,oplcolnum+2])
        for j in r_[:opl.shape[0]]:
            sid, slice1, cx1, cy1, slice2, cx2, cy2 = opl[j,:]
            target1 = ( (xyzs[:,3]==slice1) & (xyzs[:,4]==cx1) & (xyzs[:,5]==cy1) )
            if target1.any():
                subid1 = xyzs[target1,0]
            else:
                subid1 = -1
            target2 = ( (xyzs[:,3]==slice2) & (xyzs[:,4]==cx2) & (xyzs[:,5]==cy2) )
            if target2.any():
                subid2 = xyzs[target2,0]
            else:
                subid2 = -1
            tagged = vstack([tagged, r_[sid,slice1,cx1,cy1,subid1,slice2,cx2,cy2,subid2]])
        tagged = tagged[1:,:]
        savetxt(outputfile, tagged, fmt='%d\t%d\t%10.5f\t%10.5f\t%d\t%d\t%10.5f\t%10.5f\t%d')
    
def relatebacktracedopl(targetopl, refopl):
    """Relate two backtraced opl files.
    Tag the target opl_ref with the spot id of reference opl."""
    assert targetopl.endswith("_backtracedref.txt")
    assert refopl.endswith("_backtraced.txt")
    outputfile = targetopl.replace("_backtracedref","_linkoldid")    
    target = genfromtxt(targetopl)
    ref = genfromtxt(refopl)
    tagged = zeros([1,9])
    for j in r_[:target.shape[0]]:
        refid1 = refid2 = -1
        sid, slice1, cx1, cy1, subid1, slice2, cx2, cy2, subid2 = target[j,:]
        lookup1 = ((ref[:,1]==slice1) & (ref[:,4]==subid1))
        if lookup1.any():
            refid1 = ref[lookup1,0]
        else:
            refid1 = -1
        lookup2 = ((ref[:,1]==slice2) & (ref[:,4]==subid2))
        if lookup2.any():
            refid2 = ref[lookup2,0]
        else:
            refid2 = -1
        print refid1, refid2
        tagged = vstack([tagged, r_[sid,slice1,cx1,cy1,refid1,slice2,cx2,cy2,refid2]])
    tagged = tagged[1:,:]
    savetxt(outputfile, tagged, fmt='%d\t%d\t%10.5f\t%10.5f\t%d\t%d\t%10.5f\t%10.5f\t%d')        
