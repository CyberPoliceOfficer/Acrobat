import os
import PySpin
import sys
import time
import numpy as np
import cv2 as cv

# Load camera calibration
npzfile = np.load('cali_values.npz')
mtx = npzfile['mtx']
dist = npzfile['dist']
ret = npzfile['ret']
rvecs = npzfile['rvecs']
tvecs = npzfile['tvecs']

def ResizeWithAspectRatio(image, width=None, height=None, inter=cv.INTER_AREA):
    '''
        Taken from stack overflow. It's only use is for debug, it well be deleted later
    '''
    dim = None
    (h, w) = image.shape[:2]

    if width is None and height is None:
        return image
    if width is None:
        r = height / float(h)
        dim = (int(w * r), height)
    else:
        r = width / float(w)
        dim = (width, int(h * r))

    return cv.resize(image, dim, interpolation=inter)

class flir_cameras:
    def __init__(self):
        # Retrieve singleton reference to system object
        self.system = PySpin.System.GetInstance()

         # Retrieve list of cameras from the system
        self.cam_list = self.system.GetCameras()

        self.cam_ptr = None
        for i, cam in enumerate(self.cam_list):
            self.cam_ptr = cam

            # Initialize camera
            self.cam_ptr.Init()
            
            # Retrieve GenICam nodemap
            nodemap = self.cam_ptr.GetNodeMap()
            
            # In order to access the node entries, they have to be casted to a pointer type (CEnumerationPtr here)
            node_acquisition_mode = PySpin.CEnumerationPtr(nodemap.GetNode('AcquisitionMode'))

            # Retrieve entry node from enumeration node
            node_acquisition_mode_continuous = node_acquisition_mode.GetEntryByName('Continuous')

            # Retrieve integer value from entry node
            acquisition_mode_continuous = node_acquisition_mode_continuous.GetValue()
            
            # Set integer value from entry node as new value of enumeration node
            node_acquisition_mode.SetIntValue(acquisition_mode_continuous)

            #  Image acquisition must be ended when no more images are needed.
            self.cam_ptr.BeginAcquisition()
            

        
    def Acquire(self):       
        image_result = self.cam_ptr.GetNextImage(1000)

        image_converted = image_result.GetNDArray()

        image_result.Release()

        return image_converted

    def Close(self):
        # End Acquisition
        self.cam_ptr.EndAcquisition()

        # Deinitialize camera
        self.cam_ptr.DeInit()

        # The usage of del is preferred to assigning the variable to None.
        del self.cam_ptr

        # Clear camera list before releasing system
        self.cam_list.Clear()

        # Release system instance
        self.system.ReleaseInstance()

class head_pose:
    def __init__(self, head_ids = {24, 29}):
        YAH = 2
        
def main():
    camera = flir_cameras()
    img = []
    img.append(camera.Acquire())
    img.append(camera.Acquire())
    print(np.shape(img))
    camera.Close()
    cv.imshow('img', ResizeWithAspectRatio(img[1], width=1100))
    cv.waitKey()
if __name__ == '__main__':
    if main():
        sys.exit(0)
    else:
        sys.exit(1)
