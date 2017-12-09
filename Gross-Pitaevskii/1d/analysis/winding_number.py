#!/usr/bin/env python
import numpy as np

def windingNumber(f,w0=0,dPhi=np.pi):
    """
    Input:
      f  - An array of angles in radians
      w0 - Winding number of left endpoint
    
    Returns the winding number
    """
    df = np.diff(f)
    wn = np.cumsum(np.where(np.abs(df)>dPhi,-np.sign(df),0))
    return np.hstack([np.array([w0]),wn])

def unwindPhase(f):
    wn = windingNumber(f)
    return f+2.*np.pi*wn

if __name__=="__main__":
    pass
