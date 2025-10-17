# -*- coding: utf-8 -*-
import pydicom
import ctseries
class CTStudy(object):

    _meta = { 'CTDIvol':float,
              'Pitch':float,
              'Bodypart':str,
              'Software':str,
              'Single Collimation':float,
              'Total Collimation':float,
              }
    
    _tube_meta = {'Filter':str,
                  'Exposure':float,
                  'KVP':float,
                  }
    
    _type_meta = {'Kernel':str,
                  'Slice Thickness':float,
                  'RFOV':float,
                  }
    
    _tags = {'ImageType':(0x8, 0x8),
             'KVP':(0x18, 0x60),
             'RFOV':(0x18,0x1100),
             'Kernel':(0x18,0x1210),
             'study_uid':(0x20, 0xd),
             'Slice Thickness':(0x18,0x50),
             'SeriesDescription':(0x8, 0x103e),
             'Filter':(0x18, 0x1160),
             'Exposure':(0x18, 0x1152),
             'ExposureTime':(0x18, 0x1150),
             'ExposureModulation':(0x18, 0x9323),
             'Software':(0x18, 0x1020),
             'Single Collimation':(0x18,0x9306),
             'Total Collimation':(0x18,0x9307),
             'CTDIvol':(0x18, 0x9345),
             'Bodypart':(0x18, 0x15),
             'Pitch':(0x18, 0x9311),
             }

    def __init__(self, data):

        self._series = list()

        series = dict()
        sop = "1.2.840.10008.5.1.4.1.1.2"
        sop_tag = (0x8, 0x16)
        imtype_tag = (0x8, 0x8)
        seruid_tag = (0x20, 0xe)
        studuid_tag = (0x20, 0xd)

        specific_tags=[sop_tag, imtype_tag, seruid_tag, studuid_tag]

        studuid = None

        for ser in data.series_filelist:
            for im in ser:
                try:
                    dc = pydicom.dcmread(im, specific_tags=specific_tags)
                except:
                    pass
                else:
                    if dc[sop_tag].value == sop: #"CT Image Storage":
                        if studuid is None:
                            studuid = dc[studuid_tag].value
                        imtype = dc[imtype_tag].value
                        if imtype[0] == 'ORIGINAL' and imtype[1]=='PRIMARY' and imtype[2]=='AXIAL':
                            ser_uid = dc[seruid_tag].value
                            if ser_uid not in series:
                                series[ser_uid] = list()
                            series[ser_uid].append(im)
        for value in series.values():
            ser = ctseries.CTSeries(value, studuid)
            if ser.valid():
                self._series.append(ser)

    def writeStudyMetaData(self, result):
        float_res = dict()
        str_res = dict()

        for key, dtype in self._meta.items():
            tag = self._tags[key]
            try:
                value = dtype(self._series[0].getDicomValue(tag))
            except Exception as e:
                raise KeyError("Could not find dicom tag {} for study".format(tag))
                return
            if dtype == str:
                str_res[key] = value
            elif dtype == float:
                float_res[key] = value

        for series in self._series:
            tube = series.tubeSystem()
            atype = series.analysisType()
            for key, dtype in self._tube_meta.items():
                name = key + " Tube "+tube
                if ((name not in float_res) and (name not in str_res)):
                    tag = self._tags[key]
                    try:
                        value = dtype(series.getDicomValue(tag))
                    except Exception as e:
                        raise KeyError("Could not find dicom tag {} for study".format(tag))

                    if dtype == str:
                        str_res[name] = value
                    elif dtype == float:
                        float_res[name] = value
            for key, dtype in self._type_meta.items():
                name = key + " " + atype
                if ((name not in float_res) and (name not in str_res)):
                    tag = self._tags[key]
                    try:                      
                        value = dtype(series.getDicomValue(tag))
                    except Exception as e:
                        raise KeyError("Could not find dicom tag {} for study".format(tag))
                    if dtype == str:
                        str_res[name] = value.replace("[", "").replace("]","").replace("'", "").replace(",", "")
                    elif dtype == float:
                        float_res[name] = value

        for key, value in float_res.items():
            result.addFloat(key, value)
        for key, value in str_res.items():
            result.addString(key, value)

    def getSeries(self):
        return self._series
