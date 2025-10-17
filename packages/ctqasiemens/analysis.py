# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 12:27:25 2022

@author: erlean
"""

import ctstudy

import pydicom
import numpy as np
import scipy
import scipy.ndimage
import skimage.feature
import skimage.filters
import skimage.segmentation

def circle_indices2D(shape, radii, center=None, inverse=False):
    X, Y = np.indices(shape)
    if center is None:
        c = list(s//2 for s in shape)
    else:
        c=center
    if inverse:
        return (X-c[0])**2 + (Y-c[1])**2 > radii**2
    return (X-c[0])**2 + (Y-c[1])**2 <= radii**2

class AnalysisSlice(object):

    _results={'Slice thickness measured mean':float,
              'Slice thickness measured middle':float,
              'Slice thickness measured stddev':float,
              'Slice thickness comment':str,}

    def __init__(self, ctseries):
        self.results = dict()
        self.tube = ctseries.tubeSystem()
        self.calculateSliceThickness(ctseries)

    def resultsUnits(self, resname):
        if resname != 'Slice thickness comment':
            return "mm"
        return ""
    def writeResults(self, results_obj):
        for key, fn in self._results.items():
            if key in self.results:
                name = "{} Tube {} ".format(key, self.tube)
                if fn == str:
                    results_obj.addString(name, self.results[key])
                elif fn==float:
                    results_obj.addFloat(name, self.results[key])

    def calculateSliceThickness(self, ctseries):
        self.results['Slice thickness comment'] = 'Failure'

        slice_thickness = list()
        for ind, im_path in enumerate(ctseries.imagesForAnalysis()):
            dc = pydicom.dcmread(im_path)
            spacing, arr = ctseries.getArray(dc)

            ## Smooth image in direction normal to slice alu ramp
            gauss = skimage.filters.gaussian(arr, sigma=(4, 1), truncate=8.0)

            ## Finding location of ramp profile
            seed = np.argmax(gauss)
            indx, indy = np.unravel_index(seed, arr.shape)
            tolerance = gauss[indx, indy]-150
            mask = skimage.segmentation.flood(gauss, (indx, indy), tolerance=tolerance)

            ## Test for proper height of located ramp profile
            mask_height_array = np.max(mask, axis=1)
            mask_height = len(ctseries.consecutive_nonzero(mask_height_array))
            if not 9 < mask_height*spacing < 16:
                continue

            ## Center of mass for segmented ramp
            indx, indy = np.round(scipy.ndimage.center_of_mass(mask)).astype(int)

            ## Finding background signal
            background_ind = circle_indices2D(arr.shape, mask_height//2, (indx - mask_height*2, indy))
            background = arr[background_ind].mean()

            subarr = arr[indx-mask_height//4: indx+mask_height//4, indy-mask_height:indy+mask_height]

            fwhm_thres = (subarr.max() + background) * 0.5
            seed_thres = np.argmax(subarr)
            indx, indy = np.unravel_index(seed_thres, subarr.shape)

            fwhm_mask=skimage.segmentation.flood(subarr > fwhm_thres, (indx, indy), tolerance=None)
            fwhm_mask_height = len(ctseries.consecutive_nonzero(fwhm_mask.max(axis=1)))

            if not 4 < fwhm_mask_height*spacing < 16:
                continue

            fwhm = len(ctseries.consecutive_nonzero(fwhm_mask.max(axis=0)))*spacing

            slice_thickness.append((ind, fwhm * np.tan(np.deg2rad(26))))

        valid = len(slice_thickness) > 0

        if valid:
            st_arr = np.array([t for _, t in slice_thickness])
            self.results['Slice thickness measured mean'] = st_arr.mean()
            self.results['Slice thickness measured stddev'] = st_arr.std()

        if valid:
            ind_arr = np.array([ind for ind, _ in slice_thickness])
            med_ind = np.argmin(np.abs(ind_arr-np.median(ind_arr)))
            self.results['Slice thickness measured middle'] = slice_thickness[med_ind][1]

        if valid:
            self.results['Slice thickness comment'] = 'analysis successful'
        else:
            self.results['Slice thickness comment'] = 'analysis failed'


class AnalysisMTF(object):

    _results={'MTF 50 measured mean':float,
              'MTF 50 measured median':float,
              'MTF 50 measured stddev':float,
              'MTF 10 measured mean':float,
              'MTF 10 measured median':float,
              'MTF 10 measured stddev':float,
              'MTF comment':str,}

    def __init__(self, ctseries):
        self.results = dict()
        self.tube = ctseries.tubeSystem()
        self.calculateMTF(ctseries)

    def writeResults(self, results_obj):
        for key, fn in self._results.items():
            if key in self.results:
                name = "{} Tube {} ".format(key, self.tube)
                if fn == str:
                    results_obj.addString(name, self.results[key])
                elif fn==float:
                    results_obj.addFloat(name, self.results[key])

    def resultsUnits(self, resname):
        if resname != 'MTF comment':
            return "ll/mm"
        return ""
    @staticmethod
    def mtf_value(freq, mtf, value=0.5):
        for ind in range(freq.size-1):
            if mtf[ind+1] < value <= mtf[ind]:
                a = (mtf[ind+1]-mtf[ind])/(freq[ind+1]-freq[ind])
                b = mtf[ind]-a*freq[ind]
                return (value - b) / a
        return 0

    def calculateMTF(self, ctseries):
        self.results['MTF comment'] = 'Failure'

        res = list()

        sub_image_size = 2.0 # mm

        for ind, im_path in enumerate(ctseries.imagesForAnalysis()):
            dc = pydicom.dcmread(im_path)
            spacing, array = ctseries.getArray(dc)

            rad = int(sub_image_size/spacing)
            teller = 2
            while(rad < 4):
                rad = int(teller*sub_image_size/spacing)
                teller += 1

            #first wire find attempt
            peaks = skimage.feature.peak_local_max(array, min_distance=rad, threshold_abs=1000.0, num_peaks=1)
            if len(peaks) == 0:
                continue

            xc, yc = peaks[0]

            ## testing if peak is sufficent inside image
            if xc-rad*2 < 0 or xc+rad*2 >= array.shape[0]:
                continue
            if yc-rad*2 < 0 or yc+rad*2 >= array.shape[1]:
                continue

            mask = skimage.segmentation.flood(array >= array[xc, yc], (xc, yc))
            if not np.any(mask):
                continue

            cm = np.rint(scipy.ndimage.center_of_mass(array, labels=mask)).astype(int)

            x = np.squeeze(array[cm[0], cm[1]-rad:cm[1]+rad+1]).copy()
            y = np.squeeze(array[cm[0]-rad:cm[0]+rad+1, cm[1]]).copy()

            ## Detecting cupping
            cupping = False
            for line in [x, y]:
                if np.sum(line >= line.max()) > 4 and not cupping:
                    cupping = True

            fft = np.zeros(512)
            error_background = False
            for line in [x, y]:
                ## subtracting background
                background1 = line[0:rad//4].mean()
                background2 = line[-rad//4:-1].mean()
                if (background1 > -750 or background2 > -750):
                    error_background = True

                s = line - (background1+background2) / 2
                ## Normalizing PSP
                s /= s.sum()
                # zero padded FFT
                fft += np.abs(np.fft.fft(s, 512)) / 2


            if error_background:
                continue
            # FFT frequencies
            freq = np.fft.fftfreq(fft.size, spacing)
            fft = fft[freq >= 0]
            freq = freq[freq >=0]

            mtf50 = self.mtf_value(freq, fft, value=0.5)
            mtf10 = self.mtf_value(freq, fft, value=0.1)
            res.append((ind, mtf50, mtf10, cupping))


        valid = len(res) > 0

        if valid:
            mtf10  = np.array([mtf10 for _, _, mtf10, _ in res])
            mtf50  = np.array([mtf50 for _, mtf50, _, _ in res])

            self.results['MTF 50 measured mean']=np.mean(mtf50)
            self.results['MTF 10 measured mean']=np.mean(mtf10)
            self.results['MTF 50 measured median']=np.median(mtf50)
            self.results['MTF 10 measured median']=np.median(mtf10)
            self.results['MTF 50 measured stddev']=np.std(mtf50)
            self.results['MTF 10 measured stddev']=np.std(mtf10)
            cupping = any([c for _, _, _, c in res])
            if cupping:
                self.results['MTF comment'] = 'analysis successful but cupping detected'
            else:
                self.results['MTF comment'] = 'analysis successful'

        else:
            self.results['MTF comment'] = 'analysis failed'


class AnalysisHomo(object):

    _results={'Water value':float,
              'Water noise':float,
              'Uniformity 3':float,
              'Uniformity 6':float,
              'Uniformity 9':float,
              'Uniformity 12':float,
              'Uniformity mean':float,
              'Uniformity comment':str,}

    def __init__(self, ctseries):
        self.results = dict()
        self.tube = ctseries.tubeSystem()
        self.homoAnalysis(ctseries)

    def writeResults(self, results_obj):
        for key, fn in self._results.items():
            if key in self.results:
                name = "{} Tube {} ".format(key, self.tube)
                if fn == str:
                    results_obj.addString(name, self.results[key])
                elif fn==float:
                    results_obj.addFloat(name, self.results[key])
    def resultsUnits(self, resname):
        if resname != 'Uniformity comment':
            return "HU"
        return ""

    def homoAnalysis(self, ctseries):
        self.results['Uniformity comment'] = 'Failure'
        res = list()
        for ind, im_path in enumerate(ctseries.imagesForAnalysis()):
            dc = pydicom.dcmread(im_path)
            spacing, array = ctseries.getArray(dc)
            r_dist = 70.0/spacing
            radii = 10.0/spacing

            thres = array > -100.0
            cm = scipy.ndimage.center_of_mass(thres)

            water_ind = circle_indices2D(array.shape, r_dist, cm)
            water_value = np.mean(array[water_ind])
            water_noise = np.std(array[water_ind])

            uniformity_pos = {k: p for k,p  in zip([0,3,6,9,12], [(0,0), (0,1), (1,0), (0,-1), (-1,0)])}
            uniformity_vals = dict()
            test = np.copy(array)
            for k, pos in uniformity_pos.items():
                c = cm + np.array(pos)*r_dist
                ind = circle_indices2D(array.shape, radii, c)
                uniformity_vals[k] = np.mean(array[ind])
                test[ind] = test.max()

            uniformity=dict()
            for key, value in uniformity_vals.items():
                if key != 0:
                    uniformity[key] = uniformity_vals[0] - uniformity_vals[key]

            res.append((ind, water_value, water_noise, uniformity))

        valid = np.array([(-20 < v[1] < 20) for v in res])
        valid_ind = ctseries.consecutive_nonzero(valid)
        valid_res = [res[i] for i in valid_ind]

        if len(valid_res) > 0:
            i = len(valid_res) //2
            w = valid_res[i][1]
            n = valid_res[i][2]
            u = valid_res[i][3]
            self.results['Water value']=w
            self.results['Water noise']=n
            um = 0
            for key, v in u.items():
                self.results['Uniformity {}'.format(key)] = v
                um += v
            self.results['Uniformity mean'] = um / 4
            self.results['Uniformity comment'] = 'analysis successful'
        else:
            self.results['Uniformity comment'] = 'analysis failed'


def analyse(data, results, config):

    study = ctstudy.CTStudy(data)

    analysis = list()
    for ser in study.getSeries():
        if ser.valid():
            if ser.analysisType()=='SLICE':
                analysis.append(AnalysisSlice(ser))
            elif ser.analysisType()=='MTF':
                analysis.append(AnalysisMTF(ser))
            elif ser.analysisType()=='HOMO':
                analysis.append(AnalysisHomo(ser))

    study.writeStudyMetaData(results)
    # for ser in study.getSeries():
    #     if ser.valid():
    #         ser.writeSeriesMetaData(results)

    for a in analysis:
        a.writeResults(results)
