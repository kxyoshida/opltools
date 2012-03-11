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
	outputfile = "BothEnds.txt"
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

def pickupgoodtrks(oplfile, pickupfile):
    outputfile = oplfile.replace(".txt", "_gt.txt")
    opl=genfromtxt(oplfile)
    pickup=genfromtxt(pickupfile).astype('int64')
    ends = repeat(0,opl.shape[1])
    for i in r_[:pickup.size]:
        meet, = (opl[:,0]==pickup[i]).nonzero()
        if meet.size > 0:
            ends = vstack([ends, opl[meet,:]])
    ends = ends[1:,:]
    savetxt(outputfile, ends, fmt='%d\t%d\t%10.5f\t%10.5f'+'\t%d'*clip(opl.shape[1]-4,0,4))

def removeoutsidetracks(oplfile):
    """Use results table just tagged OnCell in ImageJ as an input file.
    Now adapted to Post-Track Polished opl"""
    if oplfile.find("OnCell"):
	outputfile = oplfile.replace("OnCell", "cell")
    else:
	outputfile = oplfile[:-11]+"_cell.txt"
    opl=genfromtxt(oplfile,skiprows=1)
    coln = opl.shape[1]
    idoutside=unique(opl[opl[:,-1]==0,1])
    abolish=asarray([opl[i,1] in idoutside for i in r_[:opl.shape[0]]])
    oplinside=opl[abolish==False,1:-1]
    savetxt(outputfile, oplinside, fmt='%d\t%d\t%10.5f\t%10.5f'+'\t%d'*clip(coln-6,0,4))

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

def makeptpoplreference(oplfile):
    outputfile = oplfile[:-4]+"_ref.txt"	
    opl=genfromtxt(oplfile)
    ids = unique(opl[:,0])
    ends = repeat(0,9)
    for id in ids:
	meet, = (opl[:,0]==id).nonzero()
        if meet.size > 0:
	    ends = vstack([ends, r_[opl[meet[0],:4], opl[meet[0],6], opl[meet[-1],1:4], opl[meet[-1],6]]])
    ends = ends[1:,:]
    savetxt(outputfile, ends, fmt='%d\t%d\t%10.5f\t%10.5f\t%d\t%d\t%10.5f\t%10.5f\t%d')

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

def backtrace(oplfile,polishdir="PolishedSpots"):
    outputfile = oplfile[:-4]+"_backtraced.txt"
    opl=genfromtxt(oplfile)
    xyzs=concatpolishdir(polishdir)
    oplcolnum = opl.shape[1]
    if oplcolnum == 4:
        tagged = zeros([1,oplcolnum+1])
        for j in r_[:opl.shape[0]]:
            sid, slice, cx, cy = opl[j,:]
            #            target = (xyzs[:,3]==slice)
            target = ( (xyzs[:,3]==slice) & (xyzs[:,4]==cx) & (xyzs[:,5]==cy) )
            if target.any():
                if target.sum()==1:
                    subid = xyzs[target,0]
                else:
                    subid = -target.sum()
            else:
                subid = -1
            tagged = vstack([tagged, r_[opl[j,:],subid]])
        tagged = tagged[1:,:]
        savetxt(outputfile, tagged, fmt='%d\t%d\t%10.5f\t%10.5f\t%d')

def backtraceref(oplfile,polishdir="PolishedSpots"):
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
                if target1.sum()==1:
                    subid1 = xyzs[target1,0]
                else:
                    subid1 = -target1.sum()
            else:
                subid1 = -1
            target2 = ( (xyzs[:,3]==slice2) & (xyzs[:,4]==cx2) & (xyzs[:,5]==cy2) )
            if target2.any():
                if target2.sum()==1:
                    subid2 = xyzs[target2,0]
                else:
                    subid2 = -target2.sum()
            else:
                subid2 = -1
            tagged = vstack([tagged, r_[sid,slice1,cx1,cy1,subid1,slice2,cx2,cy2,subid2]])
        tagged = tagged[1:,:]
        savetxt(outputfile, tagged, fmt='%d\t%d\t%10.5f\t%10.5f\t%d\t%d\t%10.5f\t%10.5f\t%d')
    
def relatebacktracedopl(targetopl, refopl):
    """Relate two backtraced opl files.
    Tag the target opl_ref with the spot id of reference opl."""
    assert targetopl.endswith("ref.txt")
    assert refopl.endswith("_backtraced.txt")
    outputfile = targetopl.replace("ref","linkoldid")    
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

def fixunpolished(oplfile):
    """Fill the gaps in opl track file by interpolation, which was
    produced by the failure of gaussian 2D fit in Polish_Track.java"""
    
    assert oplfile.endswith("_polish.txt")
    outputfile = oplfile.replace("_polish","_pol")
    
    opl=genfromtxt(oplfile, skiprows=1)    
    coln = opl.shape[1]
    assert coln == 8

    cids = unique(opl[:,1])
    olist = repeat(0.0,5)    

    for cid in cids:
        pind, = where(opl[:,1]==cid)
        subdata = opl[pind,:]
        subdata = subdata[subdata[:,-1]==1,:]
        if subdata.size != 0:
            subind=subdata[:,2] - subdata[0,2]
            pdata=zeros((max(subind)+1,5))
            for i in r_[:max(subind)+1]:
                if i in subind:
                    pdata[i,:]=[subdata[subind==i,j] for j in [1,2,5,6,7]]
                else:
                    headind,=nonzero([j < i for j in subind])
                    headlast=max(subind[headind])
                    tailind,=nonzero([j > i for j in subind])
                    tailfirst=min(subind[tailind])
                    gaplength=tailfirst - headlast
                    alpha=(i-headlast)*1.0/(tailfirst-headlast)
                    interpol=[cid]
                    for j in [2,5,6]:
                        interpol=r_[interpol, subdata[subind==headlast,j]*(1.0-alpha)+subdata[subind==tailfirst,j]*alpha]

                    interpol=r_[interpol,0]
                    pdata[i,:] = interpol
        olist = vstack([olist,pdata])
      
    olist = olist[1:,:]
    savetxt(outputfile, olist, fmt='%d\t%d\t%10.5f\t%10.5f\t%d')

def fix_unpolished_com(oplfile):
    """Fill the gaps in opl track file by interpolation, which was
    produced by the failure of gaussian 2D fit in Polish_Track.java"""
    
    assert oplfile.endswith("_compolish.txt")
    outputfile = oplfile.replace("_compolish","_compol")
    
    opl=genfromtxt(oplfile, skiprows=1)    
    coln = opl.shape[1]
    assert coln == 7

    cids = unique(opl[:,1])
    olist = repeat(0.0,5)    

    for cid in cids:
        pind, = where(opl[:,1]==cid)
        subdata = opl[pind,:]
        subdata = subdata[subdata[:,5]>=0,:]
        subdata = subdata[subdata[:,6]>=0,:]
        if subdata.size != 0:
            subind=subdata[:,2] - subdata[0,2]
            pdata=ones((max(subind)+1,5))
            for i in r_[:max(subind)+1]:
                if i in subind:
                    pdata[i,:-1]=[subdata[subind==i,j] for j in [1,2,5,6]]
                else:
                    headind,=nonzero([j < i for j in subind])
                    headlast=max(subind[headind])
                    tailind,=nonzero([j > i for j in subind])
                    tailfirst=min(subind[tailind])
                    gaplength=tailfirst - headlast
                    alpha=(i-headlast)*1.0/(tailfirst-headlast)
                    interpol=[cid]
                    for j in [2,5,6]:
                        interpol=r_[interpol, subdata[subind==headlast,j]*(1.0-alpha)+subdata[subind==tailfirst,j]*alpha]

                    interpol=r_[interpol,0]
                    pdata[i,:] = interpol
        olist = vstack([olist,pdata])
      
    olist = olist[1:,:]
    savetxt(outputfile, olist, fmt='%d\t%d\t%10.5f\t%10.5f\t%d')

def listpeakvelocity(oplfile):
    """Make a list of peak and mean velocities for each track
    from LifetimeTracks _gt.txt files"""
    pixelperum = 6.7
    opl=genfromtxt(oplfile)    
    
    outputfile = oplfile.replace(".txt","_peaks.txt")

    cids = unique(opl[:,0])
    olist = repeat(0,10)    

    for cid in cids:
        pind, = where(opl[:,0]==cid)
        subdata = opl[pind,:]
        #        print subdata.shape
        fri = subdata[0, 1]
        #        subind=subdata[:,2] - subdata[0,2]
        v = ((subdata[2:,2] - subdata[:-2,2])**2+(subdata[2:,3] - subdata[:-2,3])**2)**0.5
        netdist = ((subdata[-1,2] - subdata[0,2])**2+(subdata[-1,3] - subdata[0,3])**2)**0.5
        frvmax = subdata[v==max(v),1] + 1
        for frvm in frvmax:
            pdata = r_[int(oplfile[:4]), subdata[0,:4], subdata[-1,1], netdist/pixelperum, frvm, max(v)/pixelperum, mean(v)/pixelperum]
            olist = vstack([olist,pdata])

    olist = olist[1:,:]
    savetxt(outputfile, olist, fmt='%d\t%d\t%d\t%10.5f\t%10.5f\t%d\t%10.5f\t%d\t%10.5f\t%10.5f')
