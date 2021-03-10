import numpy as np

def match(xa,ya,xb,yb, threshold = 1.0):
    #pythagorean theorem, x,y is ordered pair coordinates
    #could be ra dec, though might be weird near the poles
    dx = xa - xb.reshape(xb.size,1)
    dy = ya - yb.reshape(yb.size,1)
    r  = np.sqrt(dx**2 + dy**2)
#    index_b,index_a = np.where(r < 2.0)
    index_b,index_a = np.where(r < threshold)
    return index_a,index_b


def chunk_match(xa,ya,xb,yb, threshold = 1.0):
    assert len(xa) == len(ya)
    assert len(xb) == len(yb)

    min_index = 0
    max_index = 10000
    end = xa.size

    #indices for matches
    jja = []
    jjb = []

    flag = 0
    while flag == 0:
        if max_index > end:
            sl = slice(min_index,None)
            flag =1
        else:
            sl = slice(min_index,max_index)
        xuse = xa[sl]
        yuse = ya[sl]
        
        iib,iia = match(xb,yb,xuse,yuse,
                        threshold = threshold)
        jjb = np.r_[jjb,iib ]
        jja = np.r_[jja,iia + min_index]
            
        if max_index %100000 == 0:
            print(min_index,max_index)
        min_index += 10000
        max_index += 10000
                
    jja = jja.astype(int)
    jjb = jjb.astype(int)
    dr = np.sqrt(
        (yb[jjb] - ya[jja])**2 + (xb[jjb] - xa[jja])**2
    )
                

    return jja,jjb,dr
