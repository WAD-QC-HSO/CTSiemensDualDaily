# -*- coding: utf-8 -*-
import pydicom
import numpy as np
import skimage
import scipy

class CTSeries(object):

    _types = ['MTF', 'SLICE', 'HOMO']

    _meta = { 'RFOV':float,
              'KVP':float,
              'Kernel':str,
              'Slice thickness':float,
              'Filter':str,
              'Exposure':float,
        }

    _tags = {'ImageType':(0x8, 0x8),
             'KVP':(0x18, 0x60),
             'RFOV':(0x18,0x1100),
             'Kernel':(0x18,0x1210),
             'study_uid':(0x20, 0xd),
             'Slice thickness':(0x18,0x50),
             'SeriesDescription':(0x8, 0x103e),
             'Exposure':(0x18, 0x1152),
             'Filter':(0x18, 0x1160),
             }

    def __init__(self, dc_paths : list, study_uid:str):
        self._dcs = dc_paths
        self._valid = False
        self._metaData=dict()
        self._tube="A"
        self._analysisType=""

        self.readImages(study_uid)

        self._imagesForAnalysis = list()
        self.findImagesForAnalysis()

    def analysisType(self):
       return self._analysisType
   
    def tubeSystem(self):
        return self._tube
    
    def valid(self):
        return self._valid
    
    def imagesForAnalysis(self):
        return self._imagesForAnalysis

    def readImages(self, study_uid :str):
        if len(self._dcs) == 0:
            return

        tags = [t for t in self._tags.values()]
        dc = pydicom.dcmread(self._dcs[0], specific_tags=tags)
        it = dc[(0x8, 0x8)].value
        if len(it) >= 6:
            self._tube= it[6][-1]

        if study_uid != dc[self._tags['study_uid']].value:
            return
        for key, typ in self._meta.items():
            self._metaData[key] = typ(dc[self._tags[key]].value)

        series_desc = str(dc[self._tags['SeriesDescription']].value).lower()
        if 'mtf' in series_desc or 'wire' in series_desc:
            self._analysisType='MTF'
        elif 'slice' in series_desc:
            self._analysisType='SLICE'
        elif 'homo' in series_desc or 'noise' in series_desc:
            self._analysisType='HOMO'
        else:
            self._analysisType=''

        self._valid = (self._tube in ['A', 'B']) and (self._analysisType in self._types)


        def sort_fun(x):
            dx = pydicom.dcmread(x, specific_tags=[(0x20, 0x32),])
            return dx[(0x20, 0x32)].value[2]

        self._dcs.sort(key=sort_fun)

    def findImagesForAnalysis(self):
        radiis = {'HOMO': [  98, 99, 100, 101, 102],
                  'SLICE': [78, 79, 80, 81, 82],
                  'MTF': [28, 29, 30, 31, 32],
                  }

        phantom_radi = np.array(radiis[self._analysisType])
        valid = np.zeros(len(self._dcs), dtype=int)
        for ind, impath in enumerate(self._dcs):
            dc = pydicom.dcmread(impath)
            px, arr = self.getArray(dc)

            edges = skimage.feature.canny(arr,2.0, low_threshold=50, high_threshold=400)

            pix_radi = phantom_radi / px

            hough_res = skimage.transform.hough_circle(edges, pix_radi)

            accums, cx, cy, radii = skimage.transform.hough_circle_peaks(hough_res, pix_radi, total_num_peaks=1)
            if accums.size > 0:
                valid[ind] = accums.sum() > 0.25

        valid_inds = self.consecutive_nonzero(valid)

        self._imagesForAnalysis = list([fp for ind, fp in enumerate(self._dcs) if ind in valid_inds])
        self._valid = self._valid and len(self._imagesForAnalysis) > 0

    def getDicomValue(self, tag):
        if len(self._imagesForAnalysis) > 0:
            try:
                dc = pydicom.dcmread(self._imagesForAnalysis[0], specific_tags = [tag,])
                value = dc[tag].value
            except:
                pass
            else:
                return value
        return None

    @staticmethod
    def consecutive_nonzero(array):
        ind = np.arange(array.size)
        splits = np.split(ind, ind[~array])
        split_ind = np.argmax(np.array([array[s].sum() for s in splits]))
        return np.array([i for i in ind[splits[split_ind]] if array[i]])

    @staticmethod
    def getArray( dc):
        lut = pydicom.pixel_data_handlers.apply_modality_lut
        px = sum(dc[0x28, 0x30].value)/2
        return px, lut(dc.pixel_array.astype(np.float64), dc)
