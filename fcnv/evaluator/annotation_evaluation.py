#!/usr/bin/python2

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse

normal=3

def main():
    parser = argparse.ArgumentParser(description='Evaluates FCNV prediction recall and precision. Takes an annotation file produced by FCNV.')
    parser.add_argument('filenames', type=str, nargs=1, help='path to the annotation file')
    args = parser.parse_args()

    results_file=open(args.filenames[0], "r")
    results=[]

    for line in results_file:
            line = line.strip().split(' ')
            results.append({"pos": int(line[0]),
                            "real": int(line[1]),
                            "pred": int(line[2]),})

    foundSnips=0
    regionSnips=0
    for res in results:
        if res["real"]!=normal:
            regionSnips=regionSnips+1
            if res["real"]==res["pred"]:
                foundSnips=foundSnips+1
    #print foundSnips
    #print regionSnips

    if regionSnips>0:
        recall=float(foundSnips)/regionSnips
    else:
        recall='NAN'

    
    

    foundSnips=0
    regionSnips=0
    for res in results:
       if res["pred"]!=normal:
           regionSnips=regionSnips+1
           if res["real"]==res["pred"]:
               foundSnips=foundSnips+1
               
    if regionSnips>0:
        precision=float(foundSnips)/regionSnips
    else:
        precision='NAN'
        
    #print foundSnips
    #print regionSnips
    

    print "recall: "+str(recall)
    print "precision: "+str(precision)
    #print "F-Measure: "+str(fmeasure)

    length=0
    lastCH='!'

    for res in results:
        if res["pred"]!=res["real"]:
            if length==0:
                print str(res["pred"])+": ",
            length=length+1
        if length!=0 and res["pred"]==res["real"]:
            print length
            length=0

    length=0
    lastCH='!'
    count=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    for res in results:
        if res["pred"]==lastCH:
            length=length+1
            count[res["real"]]=count[res["real"]]+1
        if res["pred"]!=lastCH:
            majority=0
            ctMax=count[0]
            i=0
            for ct in count:
                if ct>ctMax:
                    ctMax=ct
                    majority=i
                i=i+1

            if (length!=0):
                if length!=ctMax or majority!=lastCH:
                    print str(res["pos"])+' predicted: '+str(lastCH)+' ['+str(length)+"] real: "+str(majority)+" ["+str(ctMax)+"]"
            
            count=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            length=1
            lastCH=res["pred"]
            count[res["real"]]=count[res["real"]]+1

    print ''
    print 'Inside the region'
    inTheRegion=False
    for res in results:
        if inTheRegion==False and res["real"]!=normal:
            count=0
            lastch=res['pred']
            inTheRegion=True
        if inTheRegion and res["real"]==normal:
            print '[('+str(lastch)+') - '+str(count)+']',
            break
        if inTheRegion:
            if res['pred']==lastch:
                count=count+1
            else:
                print '[('+str(lastch)+') - '+str(count)+']',
                count=1
                lastch=res['pred']
            
            
            


if __name__ == '__main__':
    main()
