import numpy as np
import cv2 as cv
import glob

def ResizeWithAspectRatio(image, width=None, height=None, inter=cv.INTER_AREA):
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


videos = glob.glob('path_videos/*')
for fname in videos:
    npzfile = np.load(fname)
    images = npzfile['images']
    for img in images:
        cv.imshow(fname, ResizeWithAspectRatio(img, width=1100))
        cv.waitKey()